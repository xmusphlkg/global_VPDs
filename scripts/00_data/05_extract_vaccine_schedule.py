"""
提取疫苗免疫程序数据
记录各国各疫苗的推荐剂次和接种年龄
"""

import pandas as pd
import numpy as np

print("正在读取免疫程序数据...")
schedule_raw = pd.read_excel('data/raw/who_vaccine-schedule-data.xlsx', sheet_name='Data')

print(f"原始数据: {len(schedule_raw)}行")
print(f"列名: {schedule_raw.columns.tolist()}")

print("\n数据预览:")
print(schedule_raw.head(10))

# 查看数据结构
print("\n数据列分析:")
for col in schedule_raw.columns:
    non_null = schedule_raw[col].notna().sum()
    print(f"  {col}: {non_null}/{len(schedule_raw)} 非空")

# 保存处理后的数据
output_file = 'data/raw/vaccine_schedule.csv'
schedule_raw.to_csv(output_file, index=False)

print(f"\n✅ 免疫程序数据已保存")
print(f"   文件: {output_file}")
print(f"   总记录数: {len(schedule_raw):,}")

# 如果有国家和疫苗列
if 'ISO_3_CODE' in schedule_raw.columns or 'CODE' in schedule_raw.columns:
    country_col = 'ISO_3_CODE' if 'ISO_3_CODE' in schedule_raw.columns else 'CODE'
    print(f"   国家数: {schedule_raw[country_col].nunique()}")

# 查看疫苗相关列
vaccine_cols = [col for col in schedule_raw.columns if 'VACCINE' in col.upper() or 'ANTIGEN' in col.upper()]
if vaccine_cols:
    print(f"\n疫苗相关列: {vaccine_cols}")
    for col in vaccine_cols:
        unique_count = schedule_raw[col].nunique()
        print(f"  {col}: {unique_count}个唯一值")
