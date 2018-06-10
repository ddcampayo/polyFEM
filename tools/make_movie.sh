##  ffmpeg -r 1 -i parts%5d.png -c:v libx264 -pix_fmt yuv420p -r 2 movie.mp4
ffmpeg -r 33.2 -pattern_type glob -i "*.png"  -c:v libx264 -pix_fmt yuv420p -r 10  -vf  scale=-2:720 movie.mp4 
##ffmpeg -r 32 -pattern_type glob -i "*.png"  -vcodec libx264 -y  -an movie.mp4  -vf  scale=-2:720
