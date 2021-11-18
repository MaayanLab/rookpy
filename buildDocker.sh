
docker build -f DockerFile -t maayanlab/rookpy .
#docker push maayanlab/rookpy

docker kill rookpy
docker rm rookpy

docker run --name rookpy -e BASE_NAME="rookpy" -e TOKEN="xoxo" -e DATA="https://mssm-seq-matrix.s3.amazonaws.com/rooky_data.pkl" -p 5005:5000 -i maayanlab/rookpy
