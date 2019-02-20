# !/usr/bin/env python
# -*- coding:utf-8 -*-
import argparse
import pandas as pd
from time import clock

# 对mapseq文件处理，得到genus和species水平上的物种绝对丰度表和相对丰度表
"""
# 将所有的mapseq_merged文件转移
cp -r 04_mapseq/ERR*.trimmed.mapseq mapseq_analysis
cd mapseq_analysis
# 批量为文件改名
for name in `ls`; do echo $name ${name%.trimmed.mapseq}; done
# 批量修改文件名字
for name in `ls`; do mv $name ${name%.trimmed.mapseq}; done
# 打印输出文件名
for name in `ls`; do echo $name; done > /mnt/raid1/laisenying/毕业设计/ASD/16s/PRJEB15420/mapseq_filename.txt
# 批处理文件，删除第一行
sed -i '1d' ERR*
# 批处理文件，删除掉第一行的#
sed -i 's/#//' ERR*

执行python文件
python abundance.py -i mapseq_filename.txt \
  -g 0.4 \
  -s 0.4 \
  -o /mnt/raid1/laisenying/毕业设计/ASD/16s/PRJNA282013/result
"""

parser = argparse.ArgumentParser(description="after_mapseq_analysis")
parser.add_argument('-i', "--infile", help="the mapseq_result_filename.txt that needed analysis", default="mapseq_filename.txt")
parser.add_argument('-g', '--genus_cut', default=0.4, type=float, help="filtered thread of genus combined_cf")
parser.add_argument('-s', '--species_cut', default=0.4, type=float, help="filtered thread of species combined_cf")
parser.add_argument('-o', '--output', help="文件输出位置")
args = parser.parse_args()

t0 = clock()
# 生成绝对丰度表
f = open("%s"%args.infile, 'r')
lines = f.readlines()
file_list = []
for line in lines:
    file_list.append(line.strip())

genus_absolute_counts = []
species_absolute_counts = []
for infile in file_list:
    df = pd.read_table(infile)
    # 统计genus水平上的数据
    df_genus = df[['Genus', 'combined_cf.5']]
    print("%s 原来reads数量为："%infile + str(len(df_genus)))
    # 选择置信度大于genus_cutoff的genus
    df_genus_good = df_genus[df_genus['combined_cf.5'] > args.genus_cut]
    print("%s genus水平上过滤后reads数量为："%infile + str(len(df_genus_good)))
    genus_count = pd.value_counts(df_genus_good['Genus'])
    # 全部转换为DataFrame
    s1 = pd.Series(genus_count.index)
    s2 = pd.Series(genus_count.values)
    df_genus_counts = pd.DataFrame([s1, s2])
    df1_genus_counts = df_genus_counts.T
    df1_genus_counts.columns = ['Genus_absolute_bundance', infile]
    genus_absolute_counts.append(df1_genus_counts)
    # 统计species水平上的绝对丰度
    df_species = df[['Species', 'combined_cf.6']]
    df_species_good = df_species[df_species['combined_cf.6'] > args.species_cut]
    print("%s species水平上过滤后reads数量为："%infile + str(len(df_species_good)))
    species_count = pd.value_counts(df_species_good['Species'])
    s3 = pd.Series(species_count.index)
    s4 = pd.Series(species_count.values)
    df_species_counts = pd.DataFrame([s3, s4])
    df1_species_counts = df_species_counts.T
    df1_species_counts.columns = ['Species_absolute_bundance', infile]
    species_absolute_counts.append(df1_species_counts)


# 合并genus水平上的绝对丰度
genus_absolute_merged = genus_absolute_counts[0]
for i in range(len(genus_absolute_counts)-1):
    df_merged = pd.merge(genus_absolute_merged, genus_absolute_counts[i+1], how='outer')
    genus_absolute_merged = df_merged
# 将NA值补上
genus_absolute_merged = genus_absolute_merged.fillna(0)
genus_absolute_merged.to_csv('%s/genus_absolute_abundance.csv' % args.output)

# 合并Species水平上的绝对丰度
species_absolute_merged = species_absolute_counts[0]
for i in range(len(species_absolute_counts)-1):
    df_merged = pd.merge(species_absolute_merged, species_absolute_counts[i+1], how='outer')
    species_absolute_merged = df_merged
species_absolute_merged = species_absolute_merged.fillna(0)
species_absolute_merged.to_csv('%s/species_absolute_abundance.csv' % args.output)


#生成genues水平上的相对丰度
genus_relative_merged = genus_absolute_merged
species_relative_merged = species_absolute_merged
print("生成相对丰度表...")
for infile in file_list:
    species_relative_merged[infile] = species_relative_merged[infile]/species_relative_merged.sum()[infile]
    genus_relative_merged[infile] = genus_relative_merged[infile] / genus_relative_merged.sum()[infile]
species_relative_merged.to_csv("%s/species_relative_abundance.csv" % args.output)
genus_relative_merged.to_csv("%s/genus_relative_abundance.csv" % args.output)
t1 = clock()
print("time:"+str(t1-t0))
