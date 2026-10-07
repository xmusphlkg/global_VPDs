"""
重新提取完整的疫苗覆盖率数据
包含所有映射中需要的疫苗，优先使用WUENIC，其他数据源作为补充
"""

import pandas as pd
import numpy as np

# 读取疾病-疫苗映射
print("读取疾病-疫苗映射...")
mapping = pd.read_csv('data/processed/vaccine_disease_mapping.csv')
vaccines_needed = mapping['vaccine_code'].unique()
print(f"需要提取的疫苗数: {len(vaccines_needed)}")
print(f"疫苗列表: {sorted(vaccines_needed)}")

# 读取WHO覆盖率数据
print("\n正在读取WHO疫苗覆盖率数据...")
coverage_raw = pd.read_excel('data/raw/who_immunization_portal_coverage.xlsx', sheet_name='Data')

print(f"原始数据: {len(coverage_raw):,}行")

# 筛选需要的疫苗
coverage_filtered = coverage_raw[coverage_raw['ANTIGEN'].isin(vaccines_needed)].copy()
print(f"筛选后记录数: {len(coverage_filtered):,}")

# 优先使用WUENIC，然后OFFICIAL，最后ADMIN
# 对于PAB等特殊指标，使用其专门的category
print("\n处理数据源优先级...")
coverage_filtered['priority'] = coverage_filtered['COVERAGE_CATEGORY'].map({
    'WUENIC': 1,
    'OFFICIAL': 2,
    'ADMIN': 3,
    'PAB': 1,  # PAB有专门的category
    'HPV': 1
})

# 按国家-年份-疫苗分组，取优先级最高的数据源
coverage_filtered = coverage_filtered.sort_values('priority')
coverage_dedup = coverage_filtered.groupby(['CODE', 'YEAR', 'ANTIGEN'], as_index=False).first()

print(f"去重后记录数: {len(coverage_dedup):,}")

# 提取需要的列
df_coverage = coverage_dedup[['CODE', 'YEAR', 'ANTIGEN', 'COVERAGE', 'COVERAGE_CATEGORY']].copy()
df_coverage.columns = ['country', 'year', 'vaccine', 'coverage', 'source']

# 清理数据
df_coverage = df_coverage.dropna(subset=['country', 'year', 'vaccine'])
df_coverage['year'] = df_coverage['year'].astype(int)
df_coverage['coverage'] = pd.to_numeric(df_coverage['coverage'], errors='coerce')

# 排序
df_coverage = df_coverage.sort_values(['country', 'vaccine', 'year'])

# 保存
output_file = 'data/raw/coverage_data.csv'
df_coverage.to_csv(output_file, index=False)

print(f"\n✅ 覆盖率数据已更新")
print(f"   文件: {output_file}")
print(f"   总记录数: {len(df_coverage):,}")
print(f"   国家数: {df_coverage['country'].nunique()}")
print(f"   疫苗数: {df_coverage['vaccine'].nunique()}")
print(f"   年份范围: {df_coverage['year'].min()} - {df_coverage['year'].max()}")

print("\n各疫苗记录数:")
vaccine_counts = df_coverage.groupby('vaccine').size().sort_values(ascending=False)
for vaccine, count in vaccine_counts.items():
    print(f"  {vaccine}: {count:,}")

print("\n各数据源记录数:")
source_counts = df_coverage.groupby('source').size().sort_values(ascending=False)
for source, count in source_counts.items():
    print(f"  {source}: {count:,}")

# 检查缺失的疫苗
vaccines_in_data = set(df_coverage['vaccine'].unique())
vaccines_missing = set(vaccines_needed) - vaccines_in_data
if vaccines_missing:
    print(f"\n⚠️  映射中要求但数据中缺失的疫苗: {sorted(vaccines_missing)}")
else:
    print(f"\n✅ 所有需要的疫苗都已包含！")

print("\n数据预览:")
print(df_coverage.head(10))
