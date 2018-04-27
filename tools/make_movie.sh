##  ffmpeg -r 1 -i parts%5d.png -c:v libx264 -pix_fmt yuv420p -r 2 movie.mp4
ffmpeg -r 16 -pattern_type glob -i "part*.png"  -c:v libx264 -pix_fmt yuv420p -r 10 movie.mp4
