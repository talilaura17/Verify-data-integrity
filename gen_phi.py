import json

import os

from charm.toolbox.pairinggroup import PairingGroup, G1, ZR, pair

from functools import reduce



# 初始化雙線性配對群

group = PairingGroup('SS512')



# 將數據儲存到同一個檔案

def save_data(file_name, g, Pp, file_path, n, s, psi_values, u_values):

    data = {

        "g": group.serialize(g).hex(),

        "Pp": group.serialize(Pp).hex(),  # 將 Pp 序列化並轉換成十六進制

        "file_path": file_path,

        "n": n,

        "s": s,

        "phi_values": [group.serialize(phi).hex() for phi in psi_values],  # 將 phi 簽名轉換成十六進制字串

        "u_values": [group.serialize(u).hex() for u in u_values]  # 將 u_values 序列化

    }



    # 將數據保存到檔案

    with open(file_name, 'w') as f:

        json.dump(data, f)

    print(f"數據已保存到 {file_name}")



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



# 互動部分

def interactive_program():

    # 輸入檔案路徑

    file_path = input("請輸入檔案路徑: ")

    

    # 輸入 file_id

    file_id = input("請輸入 file_id: ")

    

    if not os.path.exists(file_path):

        print("檔案不存在，請重新執行程式並提供有效的路徑")

        return



    # 計算檔案大小和可以切成的總塊數

    file_size = os.path.getsize(file_path)

    block_size = 64  # 每個 ZR 元素為 64 字節

    total_blocks = file_size // block_size + (1 if file_size % block_size != 0 else 0)



    print(f"檔案大小為 {file_size} 字節，可以切成 {total_blocks} 個 ZR 大小的塊。")



    # 輸入想切成多少個文件塊和每個文件塊的部分數

    print(f"(hint: n*s應等於{total_blocks})")

    n = int(input("請輸入要切成的文件塊數 (n): "))

    s = int(input("請輸入每個文件塊的部分數 (s): "))



    if n * s > total_blocks:

        print("錯誤：總區塊數超過檔案可切分的區塊數，請重新執行程式並重新輸入。")

        return



    # 載入檔案並分割

    file_blocks = load_file_blocks(file_path, block_size, n, s)



    # 生成隨機私鑰 alpha

    alpha = group.random(ZR)



    # G1 群生成元

    g = group.random(G1)

    print("g:", g)



    # 公鑰 Pp = g^alpha

    Pp = g ** alpha

    print("Pp:", Pp)



    # 生成 a_j 和 u_j

    a_values = [group.random(ZR) for _ in range(s)]  # 隨機生成的 a_j 值

    u_values = [g ** a for a in a_values]  # 計算 u_j = g^a_j

    print("u_values[0]:", u_values[0])



    # 生成 phi 簽名

    phi_values = []

    for i in range(n):

        H_id_i = group.hash(file_id + str(i), G1)  # 生成哈希

        product_term = reduce(lambda x, y: x * y, [u_values[j] ** file_blocks[i][j] for j in range(s)])

        phi_i = (H_id_i * product_term) ** alpha

        phi_values.append(phi_i)

	

    print("phi_values[0]:", phi_values[0])

    # 將數據保存到 file_id.txt

    output_file = f"{file_id}.txt"

    save_data(output_file, g, Pp, file_path, n, s, phi_values, u_values)



    print("程式結束")



# 執行互動程式

interactive_program()
