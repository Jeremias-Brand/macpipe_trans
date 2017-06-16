#! /bin/bash

macpipe_dir=$1
deploy_dir=$2
# separate the things to move by space
things_to_mv="macpipe.trans.main.snake config.dev.yaml vouchers.config scripts R"
for i in $things_to_mv; do
echo cp -R  $1"/"$i $2
done


