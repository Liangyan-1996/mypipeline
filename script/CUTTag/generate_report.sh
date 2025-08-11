#!/bin/bash

# 生成CUT&RUN分析结果HTML报告

# 设置输出HTML文件
output_html="summary_report.html"

# 清空之前的HTML文件
> "$output_html"

# 写入HTML头部
cat << EOF >> "$output_html"
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CUT&RUN Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }
        h1 { color: #333; text-align: center; }
        h2 { color: #555; margin-top: 30px; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        th { background-color: #f2f2f2; font-weight: bold; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        tr:hover { background-color: #f5f5f5; }
        .image-container { display: flex; flex-wrap: wrap; gap: 20px; margin: 20px 0; }
        .image-box { flex: 1 1 300px; box-shadow: 0 0 10px rgba(0,0,0,0.1); padding: 10px; border-radius: 5px; }
        .image-box img { max-width: 100%; height: auto; cursor: pointer; }
        .modal { display: none; position: fixed; z-index: 1000; left: 0; top: 0; width: 100%; height: 100%; overflow: auto; background-color: rgba(0,0,0,0.9); }
        .modal-content { margin: auto; display: block; max-width: 90%; max-height: 90%; }
        .close { position: absolute; top: 15px; right: 35px; color: #f1f1f1; font-size: 40px; font-weight: bold; cursor: pointer; }
        .close:hover { color: #bbb; }
    </style>
</head>
<body>
    <h1>CUT&RUN Analysis Report</h1>

    <!-- 图片查看器模态框 -->
    <div id="imageModal" class="modal">
        <span class="close">&times;</span>
        <img class="modal-content" id="modalImage">
    </div>

    <h2>Quality Control Metrics</h2>
    <table>
        <tr>
            <th>Sample</th>
            <th>Raw Reads</th>
            <th>Unique Aligned</th>
            <th>Unique Aligned Rate</th>
            <th>Multi Aligned</th>
            <th>Multi Aligned Rate</th>
            <th>Overall Alignment Rate</th>
            <th>PCR Duplication Rate</th>
            <th>MT Proportion</th>
            <th>Clean Fragments</th>
            <th>FRiP</th>
        </tr>
EOF

# 遍历log目录下的所有log文件
for log_file in log/*.log; do
    # 获取样品名称
    sample=$(basename "$log_file" .log)
    echo "Processing sample: $sample"

    # 从log文件中提取QC指标
    raw_reads=$(grep -m 1 "reads; of these:" "$log_file" | awk '{print $1}')
    unique_aligned=$(grep -m 1 "aligned concordantly exactly 1 time" "$log_file" | awk '{print $1}')
    unique_aligned_rate=$(grep -m 1 "aligned concordantly exactly 1 time" "$log_file" | awk '{print $2}' | tr -d '()')
    multi_aligned=$(grep -m 1 "aligned concordantly >1 times" "$log_file" | awk '{print $1}')
    multi_aligned_rate=$(grep -m 1 "aligned concordantly >1 times" "$log_file" | awk '{print $2}' | tr -d '()')
    overall_alignment_rate=$(grep -m 1 "overall alignment rate" "$log_file" | awk '{print $1}')
    duplicates=$(grep -m 1 "found " "$log_file" | grep -o '[0-9]*' | head -1)
    end_pairs=$(grep -m 1 "sorted " "$log_file" | awk '{print $2}')
    mt_proportion=$(grep -m 1 "MT proportion:" "$log_file" | awk '{print $3}')
    clean_fragments=$(grep -m 1 "Clean Fragments:" "$log_file" | awk '{print $3}')
    frip=$(grep -m 1 "FRiP:" "$log_file" | awk '{print $2}')

    # 计算PCR重复率
    if [ -n "$duplicates" ] && [ -n "$end_pairs" ] && [ "$end_pairs" -ne 0 ]; then
        pcr_duplication_rate=$(echo "scale=5; ($duplicates / 2) * 100 / $end_pairs" | bc)
        pcr_duplication_rate="${pcr_duplication_rate}%"
    else
        pcr_duplication_rate="N/A"
    fi

    # 写入表格行
    cat << EOF >> "$output_html"
        <tr>
            <td>$sample</td>
            <td>$raw_reads</td>
            <td>$unique_aligned</td>
            <td>$unique_aligned_rate</td>
            <td>$multi_aligned</td>
            <td>$multi_aligned_rate</td>
            <td>$overall_alignment_rate%</td>
            <td>$pcr_duplication_rate</td>
            <td>$mt_proportion</td>
            <td>$clean_fragments</td>
            <td>$frip</td>
        </tr>
EOF

done

# 闭合表格
cat << EOF >> "$output_html"
    </table>

    <h2>Fragment Length Distribution</h2>
    <div class="image-container">
EOF

# 添加片段长度分布图
for length_image in mapping/*.length.png; do
    if [ -f "$length_image" ]; then
        sample=$(basename "$length_image" .length.png)
        # 读取图像文件并转换为base64
        image_data=$(base64 -w 0 "$length_image")
        cat << EOF >> "$output_html"
        <div class="image-box">
            <h3>$sample</h3>
            <img src="data:image/png;base64,$image_data" alt="$sample length distribution" onclick="openModal(this.src)">
        </div>
EOF
    fi
done

# 闭合图像容器
cat << EOF >> "$output_html"
    </div>

    <h2>Signal on Peaks</h2>
    <div class="image-container">
EOF

# 添加peak信号分布热图
for peak_image in mapping/*.png; do
    # 跳过长度分布图
    if [[ "$peak_image" != *.length.png ]] && [ -f "$peak_image" ]; then
        sample=$(basename "$peak_image" .png)
        # 读取图像文件并转换为base64
        image_data=$(base64 -w 0 "$peak_image")
        cat << EOF >> "$output_html"
        <div class="image-box">
            <h3>$sample</h3>
            <img src="data:image/png;base64,$image_data" alt="$sample peak signal" onclick="openModal(this.src)">
        </div>
EOF
    fi
done

# 闭合HTML文件
cat << EOF >> "$output_html"
    </div>

    <script>
        // 图片查看器功能
        const modal = document.getElementById("imageModal");
        const modalImg = document.getElementById("modalImage");
        const closeBtn = document.getElementsByClassName("close")[0];

        function openModal(imgSrc) {
            modal.style.display = "block";
            modalImg.src = imgSrc;
        }

        closeBtn.onclick = function() {
            modal.style.display = "none";
        }

        window.onclick = function(event) {
            if (event.target == modal) {
                modal.style.display = "none";
            }
        }
    </script>
</body>
</html>
EOF

# 使脚本可执行
chmod +x "$output_html"

# 输出完成信息
echo "HTML report generated: $output_html"