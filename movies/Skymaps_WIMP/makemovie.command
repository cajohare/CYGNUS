cd ~/Work/Zaragoza/CYGNUS/fort/movies/Skymaps_WIMP
process.StartInfo.UseShellExecute = true;
ffmpeg -r 20 -f image2 -s 1920x1080 -i %d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p -y Skymaps_WIMP.mp4
exit