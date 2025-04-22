#!/bin/bash

# 设置输出文件
output="ATAC_seq_summary.html"

# 创建HTML文件头部
cat > $output << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>ATAC-seq Analysis Summary</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: center; }
        th { background-color: #f2f2f2; }
        .section { margin: 20px 0; }
        h1, h2, h3 { color: #333; }
        .thumbnail { width: 200px; height: auto; cursor: pointer; transition: transform 0.2s; }
        .thumbnail:hover { transform: scale(1.05); }
        .modal { display: none; position: fixed; z-index: 1; left: 0; top: 0; width: 100%; height: 100%; 
                overflow: auto; background-color: rgba(0,0,0,0.9); }
        .modal-content { margin: auto; display: block; max-width: 90%; max-height: 90vh; }
        .close { color: #f1f1f1; font-size: 40px; font-weight: bold; position: absolute; right: 35px;
                top: 15px; cursor: pointer; }
    </style>
</head>
<body>
<h1>ATAC-seq Analysis Summary</h1>

<!-- 模态框 -->
<div id="imageModal" class="modal">
    <span class="close">&times;</span>
    <img class="modal-content" id="modalImg">
</div>

<script>
    var modal = document.getElementById("imageModal");
    var modalImg = document.getElementById("modalImg");
    var span = document.getElementsByClassName("close")[0];

    function showImage(imgSrc) {
        modal.style.display = "block";
        modalImg.src = imgSrc;
    }

    span.onclick = function() {
        modal.style.display = "none";
    }

    window.onclick = function(event) {
        if (event.target == modal) {
            modal.style.display = "none";
        }
    }
</script>
EOF

# 添加质控指标表格
echo "<h2>Sample Quality Control Metrics</h2>" >> $output
echo "<table>" >> $output
echo "<tr><th>Sample</th><th>MT Proportion</th><th>Clean Fragments</th><th>FRiP</th><th>TSS Enrichment Score</th></tr>" >> $output

# 遍历log文件添加质控指标
for logfile in log/*.log; do
    sample=$(basename $logfile .log)
    mt_prop=$(grep "MT proportion:" $logfile | awk '{print $3}')
    clean_frags=$(grep "Clean Fragments:" $logfile | awk '{print $3}')
    frip=$(grep "FRiP:" $logfile | awk '{print $2}')
    tss_score=$(grep -A 1 "\$TSSEscore" $logfile | tail -n 1 | awk '{print $2}')
    
    echo "<tr><td>$sample</td><td>$mt_prop</td><td>$clean_frags</td><td>$frip</td><td>$tss_score</td></tr>" >> $output
done
echo "</table>" >> $output

# 添加指标说明
cat >> $output << 'EOF'
<div class="section">
<h2>Metrics Explanation</h2>
<ul>
    <li><strong>MT Proportion</strong>: 线粒体DNA的比例，较低更好（通常<10%）</li>
    <li><strong>Clean Fragments</strong>: 最终用于分析的片段数量</li>
    <li><strong>FRiP</strong>: 落在peak区域内的reads比例，越高越好（通常>0.3）</li>
    <li><strong>TSS Enrichment Score</strong>: 转录起始位点富集得分，越高越好（通常>3）</li>
</ul>
</div>
EOF

# 创建样本名称列表
samples=()
for length_png in mapping/*.length.png; do
    sample=$(basename $length_png .length.png)
    samples+=($sample)
done

# 添加长度分布图表格
echo "<h3>Fragment Length Distribution</h3>" >> $output
echo "<table><tr>" >> $output
# 添加表头
for sample in "${samples[@]}"; do
    echo "<th>$sample</th>" >> $output
done
echo "</tr><tr>" >> $output
# 添加长度分布图
for sample in "${samples[@]}"; do
    length_png="mapping/${sample}.length.png"
    if [ -f "$length_png" ]; then
        length_base64=$(base64 -w 0 $length_png)
        echo "<td><img class='thumbnail' src='data:image/png;base64,$length_base64' onclick='showImage(this.src)' alt='Length Distribution'></td>" >> $output
    fi
done
echo "</tr></table>" >> $output

# 添加TSS热图表格
echo "<h3>TSS Enrichment Heatmap</h3>" >> $output
echo "<table><tr>" >> $output
# 添加表头
for sample in "${samples[@]}"; do
    echo "<th>$sample</th>" >> $output
done
echo "</tr><tr>" >> $output
# 添加TSS热图
for sample in "${samples[@]}"; do
    tss_png="mapping/${sample}.png"
    if [ -f "$tss_png" ]; then
        tss_base64=$(base64 -w 0 $tss_png)
        echo "<td><img class='thumbnail' src='data:image/png;base64,$tss_base64' onclick='showImage(this.src)' alt='TSS Heatmap'></td>" >> $output
    fi
done
echo "</tr></table>" >> $output

# 添加总结说明
cat >> $output << 'EOF'
<div class="section">
<h2>Analysis Summary</h2>
<ul>
    <li>片段长度分布图显示ATAC-seq特征的核小体周期性</li>
    <li>TSS热图展示开放染色质区域在转录起始位点的富集情况</li>
    <li>样品间的比较可以反映实验质量的差异和生物学意义</li>
</ul>
<p><i>点击图片可查看大图</i></p>
</div>
</body>
</html>
EOF
