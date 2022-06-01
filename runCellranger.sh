for i in N1,N2,F1,F2,Q1,Q2
do
  cellranger count --id=$i \
                 --transcriptome=/home/scBee/00ref/ameg_3.1 \
                 --fastqs=/home/scBee/01raw_data/$i \
                 --sample=$i \
                 --expect-cells=10000 \
                 --localcores=12 \
                 --localmem=32
done
