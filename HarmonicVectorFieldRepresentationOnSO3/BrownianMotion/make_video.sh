#!/usr/bin/env bash
# Erzeugt das Video aus PlotAxisAngle3d und schreibt in jedes Bild die Zeit
#   t = floor(n/30)*0.01
# wobei n die Iterationsnummer aus dem Dateinamen AxisAngle3dPic<n>.png ist.
# Die Frames werden von AnisotropicRotationalDiffusion.m direkt als PNG
# (exportgraphics, 300 dpi) geschrieben, es ist also keine EPS-Konvertierung
# mehr notwendig.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIR="$ROOT_DIR/PlotAxisAngle3d"
OUT_FILE="$ROOT_DIR/RotationalDiffusion.mp4"
FPS=20
FONT=DejaVu-Serif
OUT_WIDTH=1400      # Breite des Videos; die Frames sind deutlich groesser
POINTSIZE_REL=0.039 # Schriftgroesse relativ zur Bildbreite
OFFSET_X_REL=0.064  # Position der Zeitangabe (unten rechts), relativ zur Breite
OFFSET_Y_REL=0.050

mapfile -t PNG_FILES < <(find "$DIR" -maxdepth 1 -type f -name "AxisAngle3dPic*.png" | sort -V)

if [[ ${#PNG_FILES[@]} -eq 0 ]]; then
    echo "Keine PNG-Dateien in $DIR gefunden." >&2
    exit 1
fi

echo "Gefundene PNG-Dateien: ${#PNG_FILES[@]}"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

# Schriftgroesse, Offset und Zielgroesse aus dem ersten Frame bestimmen.
# Die exportgraphics-Frames sind nicht alle exakt gleich gross, daher wird
# jeder Frame auf dieselbe Leinwand OUT_WIDTH x OUT_HEIGHT gebracht.
read -r WIDTH HEIGHT < <(identify -format '%w %h\n' "${PNG_FILES[0]}")
POINTSIZE=$(awk -v w="$WIDTH" -v r="$POINTSIZE_REL" 'BEGIN{printf "%d", w*r}')
OFFSET=$(awk -v w="$WIDTH" -v x="$OFFSET_X_REL" -v y="$OFFSET_Y_REL" \
             'BEGIN{printf "+%d+%d", w*x, w*y}')
# Hoehe passend zum Seitenverhaeltnis, auf gerade Zahl gerundet (libx264/yuv420p)
OUT_HEIGHT=$(awk -v w="$WIDTH" -v h="$HEIGHT" -v ow="$OUT_WIDTH" \
                 'BEGIN{printf "%d", int(ow*h/w/2)*2}')
GEOM="${OUT_WIDTH}x${OUT_HEIGHT}"

echo "Videogroesse: $GEOM, Schriftgroesse: $POINTSIZE, Offset: $OFFSET"

j=1
for PNG_SRC in "${PNG_FILES[@]}"; do
    BASE="$(basename "$PNG_SRC" .png)"
    N="${BASE#AxisAngle3dPic}"

    # t = floor(n/30) * 0.01  (ganzzahlig gerechnet, um Rundungsfehler zu vermeiden)
    HUNDREDTHS=$(( N / 30 ))
    LABEL=$(printf 't = %d.%02d' $(( HUNDREDTHS / 100 )) $(( HUNDREDTHS % 100 )))

    PNG_FILE="$TMP_DIR/frame$(printf '%05d' "$j").png"

    echo "  Frame $j: $BASE  ->  $LABEL"

    convert "$PNG_SRC" \
       -background white -alpha remove -alpha off \
       -font "$FONT" -pointsize "$POINTSIZE" -fill black \
       -gravity SouthEast -annotate "$OFFSET" "$LABEL" \
       -resize "$GEOM" -gravity center -extent "$GEOM" \
       "$PNG_FILE"

    j=$((j+1))
done

echo "Erzeuge Video: $OUT_FILE"

ffmpeg -y -framerate "$FPS" \
  -i "$TMP_DIR/frame%05d.png" \
  -c:v libx264 \
  -pix_fmt yuv420p \
  "$OUT_FILE"

echo "Fertig: $OUT_FILE"
