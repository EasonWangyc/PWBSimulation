"""
analyze_results.py — 仿真结果后处理
读取已有的扫描结果文件，绘制 loss 热力图。
无需启动 FDTD，可独立运行。
"""

from pwb_core import plot_Ttotal_loss_heatmap_rR

# r/R 扫描结果热力图
data_file_rR = "D:/simulation/Simulation Project/simulation/PD-PWB-SMF/results/r_R_scan/T_total_sweep_results.txt"
loss_matrix, r_vals, R_vals = plot_Ttotal_loss_heatmap_rR(data_file_rR)
print("r values (μm):", r_vals)
print("R values (μm):", R_vals)

# 如需分析 h/R 扫描结果，取消注释以下代码
# data_file_hR = "D:/simulation/Simulation Project/simulation/PD-PWB-SMF/results/h_R_scan/T_total_sweep_results.txt"
# loss_matrix_hR, R_vals_hR, h_vals_hR = plot_Ttotal_loss_heatmap_rR(data_file_hR)
