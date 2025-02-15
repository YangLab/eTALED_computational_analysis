mkdir 4-0_radar

parallel -j 4 --rpl ",, s/.*\/(.*).R1.30M.fq.gz/\1/" '/data/bioind1/xingma/Projects/BEhuman/RNA_radar/20231205/RADAR_dev/RADAR -1 {} -2 {= s/.R1./.R2./ =} -o ./ -n ,, -t 20  -g hg38 -c /data/bioind1/xingma/Projects/BEhuman/RNA_radar/20231205/RADAR_dev/src/RADAR.conf &> ,,.log' ::: ../3-0_seqtk/*.R1.30M.fq.gz



