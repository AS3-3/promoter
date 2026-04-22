import pandas as pd
from collections import defaultdict

# 输入和输出文件名
input_file = r'gene-location-elements.txt'
middle_file = r'0413-40.xlsx'
output_file = r'0413-40-矩阵.xlsx'

# 设定的距离可以在这里更改，不更改则为前后60bp
distance = 40

# 读取数据文件
df = pd.read_csv(input_file, sep='\t', header=None, names=['基因id', '位置', '功能元件'])

# 获取所有不重复的基因ID
all_genes = df['基因id'].unique().tolist()

# 按功能元件分组，记录每行的原始索引
grouped = defaultdict(list)
for idx, row in df.iterrows():
    grouped[row['功能元件']].append((row['基因id'], row['位置'], idx))

output = []
global_rank = 1  # 全局名次计数

# 核心处理函数（修改为蔓延聚类逻辑）
def process_step(elements):
    global global_rank
    sub_output = []

    for func, items in elements.items():
        # 如果该功能元件下没有行，跳过
        if not items:
            continue

        # 按位置排序，便于线性扫描聚类
        items_sorted = sorted(items, key=lambda x: x[1])  # x[1] 是位置

        # 初始化第一个簇
        clusters = []
        current_cluster = [items_sorted[0]]

        # 线性扫描，根据距离判断是否属于同一簇
        for i in range(1, len(items_sorted)):
            prev_item = current_cluster[-1]
            curr_item = items_sorted[i]
            # 如果当前行位置与前一个簇尾行位置差 <= distance，则加入当前簇
            if curr_item[1] - prev_item[1] <= distance:
                current_cluster.append(curr_item)
            else:
                # 否则，保存当前簇，并开始新簇
                clusters.append(current_cluster)
                current_cluster = [curr_item]
        # 添加最后一个簇
        clusters.append(current_cluster)

        # 为每个簇分配名次并生成记录
        for cluster in clusters:
            rank = global_rank
            global_rank += 1

            # 记录簇内出现的基因（标记为1）
            cluster_genes = set()
            for gene_id, pos, idx in cluster:
                sub_output.append((gene_id, func, pos, rank, idx, 1))
                cluster_genes.add(gene_id)

            # 对于该功能元件下未在此簇出现的基因，补充类型0记录
            # 注意：这里我们仍然以全局基因列表 all_genes 来补0，
            # 保持与原始代码行为一致（如需仅对参与该元件的基因补0可调整）
            missing_genes = [g for g in all_genes if g not in cluster_genes]
            for gene_id in missing_genes:
                sub_output.append((gene_id, func, None, rank, -1, 0))

    output.extend(sub_output)

# 执行聚类处理
process_step(grouped)

# 输出第一张详细表
output_df = pd.DataFrame(output, columns=['基因id', '功能元件', '位置', '名次', '行索引', '类型'])
output_df.to_excel(middle_file, index=False, engine='openpyxl')
print("处理完成，结果已保存到", middle_file)

# 读取数据并生成矩阵
df = pd.read_excel(middle_file)

# 按名次和基因id去重，取第一行
unique_df = df.groupby(['名次', '基因id']).first().reset_index()

# 创建基因 × 名次的矩阵
matrix = unique_df.pivot(index='基因id', columns='名次', values='类型').fillna(0).astype(int)

# 保存矩阵
matrix.to_excel(output_file)
print("全部处理完成！矩阵已保存到", output_file)