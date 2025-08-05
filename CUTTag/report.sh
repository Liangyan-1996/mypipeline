#!/bin/bash

# --- 配置 ---
# 定义源目录和输出目录
OUTPUT_DIR="log"
MAPPING_DIR="mapping"
MOTIF_DIR="motif"
REPORT_HTML="$OUTPUT_DIR/summary_report.html"

# --- 检查 ---
# 检查必要的目录是否存在
if [ ! -d "$OUTPUT_DIR" ] || [ ! -d "$MAPPING_DIR" ] || [ ! -d "$MOTIF_DIR" ]; then
    echo "错误: 脚本必须在包含 'log', 'mapping', 和 'motif' 子目录的项目根目录下运行。"
    exit 1
fi

# --- 初始化 ---
# 清理旧的报告和所有软链接，以防重复运行
# 'find -type l -delete' 会删除目录下的所有软链接
rm -f "$REPORT_HTML"
find "$OUTPUT_DIR" -maxdepth 1 -type l -delete
echo "旧报告和软链接已清理。"

# --- 开始生成HTML文件 ---
# 写入HTML头部和CSS样式
cat <<EOF > "$REPORT_HTML"
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <title>CUT&RUN 分析报告</title>
    <style>
        body { font-family: sans-serif; margin: 2em; }
        h1, h2 { color: #333; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        a { color: #0066cc; text-decoration: none; }
        a:hover { text-decoration: underline; }
        .image-gallery { display: flex; flex-wrap: wrap; gap: 10px; margin-top: 10px; }
        .image-gallery img { max-height: 200px; border: 1px solid #ccc; cursor: pointer; transition: transform 0.2s; }
        .image-gallery img:hover { transform: scale(1.05); }
    </style>
</head>
<body>
    <h1>CUT&RUN 分析结果汇总</h1>

    <h2>核心指标与Motif富集分析</h2>
    <table>
        <thead>
            <tr>
                <th>样品名称</th>
                <th>唯一比对数</th>
                <th>唯一比对率</th>
                <th>PCR 重复率</th>
                <th>MT 比例</th>
                <th>Clean Fragments</th>
                <th>FRiP Score</th>
                <th>De Novo Motif</th>
                <th>Known Motif</th>
            </tr>
        </thead>
        <tbody>
EOF

# --- 核心处理逻辑：循环读取log文件 ---
# 寻找所有的log文件并逐一处理
for log_file in "$OUTPUT_DIR"/*.log; do
    # 从文件名提取样品名
    sample_name=$(basename "$log_file" .log)
    echo "正在处理样品: $sample_name"

    # --- 1. 从 log 文件提取指标 ---
    unique_line=$(grep "aligned concordantly exactly 1 time" "$log_file")
    unique_count=$(echo "$unique_line" | awk '{print $1}')
    unique_rate=$(echo "$unique_line" | awk '{print $2}' | sed 's/[()]//g')

    duplicates=$(grep "found" "$log_file" | grep "duplicates" | awk '{print $2}')
    end_pairs=$(grep "sorted" "$log_file" | grep "end pairs" | awk '{print $2}')
    pcr_dup_rate="N/A"
    if [[ -n "$duplicates" && -n "$end_pairs" && "$end_pairs" -ne 0 ]]; then
        pcr_dup_rate=$(echo "scale=4; (($duplicates / 2) / $end_pairs) * 100" | bc | awk '{printf "%.2f%%", $1}')
    fi

    # 修正：直接获取MT比例，不再乘以100
    mt_proportion=$(grep "MT proportion" "$log_file" | awk '{print $3}')

    clean_fragments=$(grep "Clean Fragments" "$log_file" | awk '{print $3}')
    frip_score=$(grep "FRiP" "$log_file" | awk '{print $2}')

    # --- 2. 处理 Motif 结果 ---
    sample_motif_dir="$MOTIF_DIR/$sample_name"
    denovo_link="N/A"
    known_link="N/A"
    
    if [ -d "$sample_motif_dir" ]; then
        echo "  > 发现 Motif 结果目录，创建软链接..."
        # 在 log 目录下，创建指向 ../motif/SAMPLE 的软链接
        ln -s "../$sample_motif_dir" "$OUTPUT_DIR/$sample_name"

        # 生成HTML超链接，路径相对于 summary_report.html
        if [ -f "$sample_motif_dir/homerResults.html" ]; then
            denovo_link="<a href=\"$sample_name/homerResults.html\" target=\"_blank\">查看报告</a>"
        fi
        if [ -f "$sample_motif_dir/knownResults.html" ]; then
            known_link="<a href=\"$sample_name/knownResults.html\" target=\"_blank\">查看报告</a>"
        fi
    fi

    # --- 3. 将所有数据写入HTML表格的一行 ---
    cat <<EOF >> "$REPORT_HTML"
            <tr>
                <td>$sample_name</td>
                <td>${unique_count:-N/A}</td>
                <td>${unique_rate:-N/A}</td>
                <td>$pcr_dup_rate</td>
                <td>${mt_proportion:-N/A}</td>
                <td>${clean_fragments:-N/A}</td>
                <td>${frip_score:-N/A}</td>
                <td>$denovo_link</td>
                <td>$known_link</td>
            </tr>
EOF
done

# --- 结束表格 ---
cat <<EOF >> "$REPORT_HTML"
        </tbody>
    </table>
EOF

# --- 处理图片 (此部分逻辑不变) ---
# 整理片段长度分布图 (*.length.png)
echo "正在整理片段长度分布图..."
cat <<EOF >> "$REPORT_HTML"
    <h2>片段长度分布</h2>
    <div class="image-gallery">
EOF
for img_path in "$MAPPING_DIR"/*_length.png; do
    if [ -f "$img_path" ]; then
        img_name=$(basename "$img_path")
        ln -s "../$MAPPING_DIR/$img_name" "$OUTPUT_DIR/$img_name"
        echo "        <a href=\"$img_name\" target=\"_blank\"><img src=\"$img_name\" alt=\"$img_name\"></a>" >> "$REPORT_HTML"
    fi
done
cat <<EOF >> "$REPORT_HTML"
    </div>
EOF

# 整理Peak信号分布热图 (*.png, 但排除 *_length.png)
echo "正在整理Peak信号热图..."
cat <<EOF >> "$REPORT_HTML"
    <h2>Peak信号热图</h2>
    <div class="image-gallery">
EOF
for img_path in $(ls "$MAPPING_DIR"/*.png 2>/dev/null | grep -v '_length.png'); do
    if [ -f "$img_path" ]; then
        img_name=$(basename "$img_path")
        ln -s "../$MAPPING_DIR/$img_name" "$OUTPUT_DIR/$img_name"
        echo "        <a href=\"$img_name\" target=\"_blank\"><img src=\"$img_name\" alt=\"$img_name\"></a>" >> "$REPORT_HTML"
    fi
done
cat <<EOF >> "$REPORT_HTML"
    </div>
EOF

# --- 结束HTML文件 ---
cat <<EOF >> "$REPORT_HTML"
</body>
</html>
EOF

echo ""
echo "报告生成完毕！"
echo "HTML报告路径: $(pwd)/$REPORT_HTML"
echo "所有依赖的图片和Motif结果目录已通过软链接方式放入 '$OUTPUT_DIR' 目录。"