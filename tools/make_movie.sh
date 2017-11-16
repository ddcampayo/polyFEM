##  ffmpeg -r 1 -i parts%5d.png -c:v libx264 -pix_fmt yuv420p -r 2 movie.mp4
ffmpeg -r 1 -pattern_type glob -i "*.png"  -c:v libx264 -pix_fmt yuv420p -r 4 movie.mp4
