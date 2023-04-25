#!/bin/bash

## install aws-cli (conda-forge)
## install tower cli (https://help.tower.nf/22.4/cli/#availability)
    # export TOWER_ACCESS_TOKEN=<access-token> --> added to .zshrc so no need to do this again


echo "Run ID: $1"

# tw -o json runs view -w <workspace_id> -i <run_id> task -t <task_id> --resources-requested | jq ".resources.machineType"
# default workspace used, so removed -w

# BWA 2.7.2: 4Sr1TIAwZ3nisp

# page_counter=1
# page_content=""
# task_ids=()
# while [ "$page_content" != "[ ]" ]
# do
#     page_content=$(tw -o json runs view -i $1 tasks --page $page_counter)
#     task_ids+=$(tw -o json runs view -i $1 tasks --page $page_counter | jq '.[].taskId')
#     task_ids+=" "
#     page_counter=$((page_counter+1))
# done

# for task_id in $task_ids
# do
#     echo "$task_id"
#     echo "$(tw -o json runs view -i $1 task -t $task_id --resources-requested --execution-time | jq -r "[.resources.machineType,.resources.cloudZone,.times.start,.times.end]  | @tsv")" >> resource.tsv
# done

exec 3< "resource.tsv"
# Loop over the rows in the file
while read -r machineType cloudZone start end; do


    # Do something with the values in the row
    echo "Column 1: $machineType"
    echo "Column 2: $cloudZone"
    echo "Column 3: $start"
done <&3

# Close the file
exec 3<&-

