avconv -r 20 -i parts%1d.png -c:v libx264 -pix_fmt yuv420p   movie.mp4
