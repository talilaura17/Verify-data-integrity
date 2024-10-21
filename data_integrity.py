from charm.toolbox.pairinggroup import PairingGroup, ZR, G1, pair

from functools import reduce



# 初始化一個雙線性配對群

group = PairingGroup('SS512')



# 生成隨機私鑰 alpha

alpha = group.random(ZR)



# G1 群生成元

g = group.random(G1)



# 公鑰 Pp = g^alpha

Pp = g ** alpha



# 假設我們有 n 個文件塊，每個文件塊有 s 個部分 (sectors)

n, s = 5, 3  # n 是文件塊數量，s 是每個塊的部分數量

file_id = "file123"  # 假設這是文件的唯一標識符



# 生成文件塊和部分數據（隨機模擬數據）

file_blocks = [[group.random(ZR) for _ in range(s)] for _ in range(n)]  # 模擬每個文件塊內的數據部分



# 生成隨機的 a_j 並計算 u_j = g^a_j

a_values = [group.random(ZR) for _ in range(s)]  # 隨機生成的 a_j 值

u_values = [g ** a for a in a_values]  # 計算 u_j = g^a_j



# 對每個文件塊生成簽名 psi_i

psi_values = []

for i in range(n):

    # 生成簽名 psi_i = (H(id_i) * 乘積(u_j^f_ij))^alpha

    H_id_i = group.hash(file_id + str(i), G1)  # 生成文件塊唯一的哈希值

    product_term = reduce(lambda x, y: x * y, [u_values[j] ** file_blocks[i][j] for j in range(s)])

    psi_i = (H_id_i * product_term) ** alpha

    psi_values.append(psi_i)



# 假設挑戰包括一些隨機選擇的文件塊

chal = [0, 2]  # 挑戰的文件塊索引

c_values = [group.random(ZR) for _ in chal]  # 挑戰的權重系數



# 聚合簽名 ψ 和線性組合 χ 的計算

psi_agg = group.init(G1, 1)  # 初始化聚合簽名

chi_agg = [0 for _ in range(s)]  # 用來存儲每個部分的線性組合結果



for i in range(len(chal)):

    psi_agg *= psi_values[chal[i]] ** c_values[i]  # 聚合簽名 ψ

    for j in range(s):

        chi_agg[j] += file_blocks[chal[i]][j] * c_values[i]  # 計算每個部分的線性組合 χ



# 驗證步驟

# 左側配對計算 lhs = e(ψ_agg, g)

lhs = pair(psi_agg, g)



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

