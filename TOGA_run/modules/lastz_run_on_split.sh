#! /bin/bash

target=$1
query=$2
target_region = $3
query_region = $4
target_info=$4
query_info=$5
output=$6

mkdir -p "${output}/tmp_lastz"

parallel -j 4 --xapply "lastz \"${target}/{1}\" \"{query}/{2}\" --format=axt | axtChain -linearGap=loose -minScore=1000 stdin ${target} ${query} ${ouput}/tmp_lastz/{}.chain" ::: "A B C" ::: "1 2 3"

# 并行运行lastz和axtChain命令
cut -f1 "${query_info}" | parallel -j 4 "lastz \"${target}/${target_region}\" \"${query}/{}\" --format=axt | axtChain -linearGap=loose -minScore=1000 stdin ${target} ${query} ${output}/tmp_lastz/{}.chain"

# 合并所有chain文件为一个文件
find "${output}/tmp_lastz/" -name '*.chain' | chainMergeSort -inputList=stdin > "${output}/lastz.chain"

# 删除临时文件夹
rm -rf "${output}/tmp_lastz"