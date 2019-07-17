cd ~/Work/Zaragoza/CYGNUS/fort/movies/Skymaps_SN
process.StartInfo.UseShellExecute = true;
ffmpeg -r 20 -f image2 -s 1920x1080 -i 1_%d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p -y -vf 'scale=-2:min(1080\,trunc(ih/2)*2)' Skymaps_SN1.mp4 
exit