"""
提取疫苗引入数据（修正版）
WHO数据格式: 国家-年份-疫苗描述-是否引入
"""

import pandas as pd
import numpy as np

print("正在读取疫苗引入数据...")
intro_raw = pd.read_excel('data/raw/who_vaccine-introduction-data.xlsx', sheet_name='Data')

print(f"原始数据: {len(intro_raw)}行")

# 筛选已引入的记录
intro_yes = intro_raw[intro_raw['INTRO'] == 'Yes'].copy()
print(f"已引入记录: {len(intro_yes)}行")

# 提取需要的列
df_intro = intro_yes[['ISO_3_CODE', 'DESCRIPTION', 'YEAR']].copy()
df_intro.columns = ['country', 'vaccine_description', 'intro_year']

# 清理数据
df_intro = df_intro.dropna(subset=['country', 'vaccine_description', 'intro_year'])
df_intro['intro_year'] = pd.to_numeric(df_intro['intro_year'], errors='coerce')
df_intro = df_intro.dropna(subset=['intro_year'])

# 每个国家-疫苗组合保留最早引入年份
df_intro = df_intro.sort_values('intro_year').groupby(['country', 'vaccine_description'], as_index=False).first()

# 创建简化的疫苗代码映射
vaccine_mapping = {
    'Measles-containing vaccine 2nd dose': 'MCV2',
    'Rubella vaccine': 'RCV1',
    'Mumps vaccine': 'Mumps',
    'Japanese Encephalitis': 'JAPENC',
    'YF (Yellow fever) vaccine': 'YFV',
    'Meningococcal meningitis vaccines (any strain)': 'Meningococcal',
    'Typhoid conjugate vaccine': 'TYPHOID_CONJ',
    'IPV (Inactivated polio vaccine)': 'IPV',
    'IPV (Inactivated polio vaccine) 2nd dose': 'IPV2',
    'Hepatitis B vaccine': 'HepB',
    'Hib (Haemophilus influenzae type B) vaccine': 'Hib',
    'HPV (Human papillomavirus) vaccine': 'HPV',
    'PCV (Pneumococcal conjugate vaccine)': 'PCV',
    'Rotavirus vaccine': 'Rotavirus',
}

df_intro['vaccine_code'] = df_intro['vaccine_description'].map(vaccine_mapping)

# 保存完整数据
output_file = 'data/raw/vaccine_introduction.csv'
df_intro.to_csv(output_file, index=False)

print(f"\n✅ 疫苗引入数据已提取")
print(f"   文件: {output_file}")
print(f"   总记录数: {len(df_intro):,}")
print(f"   国家数: {df_intro['country'].nunique()}")
print(f"   疫苗类型数: {df_intro['vaccine_description'].nunique()}")
print(f"   引入年份范围: {int(df_intro['intro_year'].min())} - {int(df_intro['intro_year'].max())}")

print("\n各疫苗引入记录数:")
vaccine_counts = df_intro.groupby('vaccine_description').size().sort_values(ascending=False)
for vaccine, count in vaccine_counts.items():
    print(f"  {vaccine}: {count}")

print("\n关键疾病的引入数据覆盖:")
key_vaccines = ['Measles-containing vaccine 2nd dose', 'Rubella vaccine', 'Mumps vaccine',
                'Japanese Encephalitis', 'YF (Yellow fever) vaccine']
for vac in key_vaccines:
    count = len(df_intro[df_intro['vaccine_description'] == vac])
    print(f"  {vac}: {count}个国家")

print("\n数据预览:")
print(df_intro.head(10))
