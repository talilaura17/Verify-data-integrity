import json

import os

from charm.toolbox.pairinggroup import PairingGroup, G1, ZR, pair

from functools import reduce



# 初始化雙線性配對群

group = PairingGroup('SS512')


# 從檔案讀取數據

def load_data(file_name):

    with open(file_name, 'r') as f:

        data = json.load(f)


    # 反序列化 g, Pp, phi_values 和 u_values

    g = group.deserialize(bytes.fromhex(data["g"]))

    Pp = group.deserialize(bytes.fromhex(data["Pp"]))

    phi_values = [group.deserialize(bytes.fromhex(phi)) for phi in data["phi_values"]]

    u_values = [group.deserialize(bytes.fromhex(u)) for u in data["u_values"]]

    file_path = data["file_path"]

    n = data["n"]

    s = data["s"]


    return g, Pp, file_path, n, s, phi_values, u_values



# 載入檔案並分割成文件塊

def load_file_blocks(file_path, block_size, n, s):

    file_blocks = []

    with open(file_path, 'rb') as file:

        for i in range(n):

            block_row = []

            for j in range(s):

                block_data = file.read(block_size)

                if len(block_data) < block_size:  # 如果不足 64 字節，則補齊

                    block_data = block_data.ljust(block_size, b'\0')

                # 將檔案塊數據轉換成 ZR 元素

                block_zr_value = group.hash(block_data, ZR)

                block_row.append(block_zr_value)

            file_blocks.append(block_row)

    return file_blocks





# 調用 load_data 函數來讀取存儲的數據

file_id = input("請輸入要讀取的 file_id: ")

file_name = f"{file_id}.txt"

g, Pp, file_path, n, s, phi_values, u_values = load_data(file_name)

block_size = 64

file_blocks = load_file_blocks(file_path, block_size, n, s)


# 挑戰包括一些隨機選擇的文件塊
print("n=",n)

# 挑戰的文件塊索引

chal_input = input("請輸入挑戰的文件塊索引（用逗號分隔，如: 0, 2）(<n): ")
chal = list(map(int, chal_input.split(',')))

c_values = [group.random(ZR) for _ in chal]  # 挑戰的權重系數

print("隨機挑戰:", c_values)



# 聚合簽名 ψ 和線性組合 χ 的計算

phi_agg = group.init(G1, 1)  # 初始化聚合簽名
chi_agg = [0 for _ in range(s)]  # 用來存儲每個部分的線性組合結果

for i in range(len(chal)):

    phi_agg *= phi_values[chal[i]] ** c_values[i]  # 聚合簽名 ψ

    for j in range(s):

        chi_agg[j] += file_blocks[chal[i]][j] * c_values[i]  # 計算每個部分的線性組合 χ







# 驗證步驟

# 左側配對計算 lhs = e(ψ_agg, g)

lhs = pair(phi_agg, g)



# 右側配對計算 rhs = e( H(id_i)^c_i , Pp ) * 乘積( e( u_j^chi_agg_j , Pp ) )

rhs = pair(group.init(G1, 1), Pp)  # 初始化右側配對計算


for i in range(len(chal)):

    H_id_i = group.hash(file_id + str(chal[i]), G1)  # 針對每個挑戰計算 H(id_i)

    rhs *= pair(H_id_i ** c_values[i], Pp)



# 計算 e(u_j^chi_agg_j, Pp)

for j in range(s):

    rhs *= pair(u_values[j] ** chi_agg[j], Pp)


# 檢查配對結果是否相同

print("左側配對結果:", lhs)

print("右側配對結果:", rhs)

print("驗證通過:", lhs == rhs)

