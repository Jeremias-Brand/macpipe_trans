#! /bin/bash
if [[ $# -ne 2 ]]; then
    echo -e "USAGE:\n deploy.sh macpipe_dir target_dir"
    exit 1
fi

macpipe_dir=$1
deploy_dir=$2
# separate the things to move by space
things_to_mv="macpipe.snakefile config.yaml vouchers.config scripts R dbs"
# echo the commands that will be executed
for i in $things_to_mv; do
echo cp -R  $1$i $2
done
# check if the user is sure about command and run
read -r -p "Deploy files? [y/N] " response
if [[ "$response" =~ ^([yY][eE][sS]|[yY])+$ ]] # allow yes and single letter
then
    for i in $things_to_mv; do
        cp -R  $1$i $2
    done
else
    echo "Aborted!"
fi


