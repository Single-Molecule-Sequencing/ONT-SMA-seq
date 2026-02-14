/**
 * D3.js ROC curve and threshold recommendation visualization.
 *
 * Provides:
 *   renderROC(containerId, data, options) - ROC curve plot
 *   renderThresholdTable(containerId, recommendations) - metrics table
 *
 * ROC data format:
 *   { fpr: [...], tpr: [...], thresholds: [...], auc: float }
 *
 * Options (ROC):
 *   { title, recommendedThreshold, onPointClick(threshold) }
 *
 * Recommendations format:
 *   [ { column, threshold, sensitivity, specificity, f1, fpr,
 *       true_positives, false_positives, reads_affected } ]
 */


/**
 * Render a ROC curve with shaded AUC area, diagonal reference, and
 * a recommended threshold marker.
 */
function renderROC(containerId, data, options) {
  const container = document.getElementById(containerId);
  if (!container) return;
  container.innerHTML = '';

  options = options || {};

  const fpr = data.fpr || [];
  const tpr = data.tpr || [];
  const thresholds = data.thresholds || [];
  const auc = data.auc;

  if (fpr.length === 0) {
    container.innerHTML = '<p>No ROC data available.</p>';
    return;
  }

  const margin = { top: 30, right: 30, bottom: 50, left: 60 };
  const width = Math.max(container.clientWidth || 500, 400) - margin.left - margin.right;
  const height = 400 - margin.top - margin.bottom;

  const xScale = d3.scaleLinear().domain([0, 1]).range([0, width]);
  const yScale = d3.scaleLinear().domain([0, 1]).range([height, 0]);

  const svg = d3.select(container)
    .append('svg')
    .attr('width', width + margin.left + margin.right)
    .attr('height', height + margin.top + margin.bottom);

  const g = svg.append('g')
    .attr('transform', `translate(${margin.left},${margin.top})`);

  // Axes
  g.append('g')
    .attr('class', 'kde-axis')
    .attr('transform', `translate(0,${height})`)
    .call(d3.axisBottom(xScale).ticks(10));

  g.append('g')
    .attr('class', 'kde-axis')
    .call(d3.axisLeft(yScale).ticks(10));

  // Axis labels
  g.append('text')
    .attr('x', width / 2)
    .attr('y', height + 40)
    .attr('text-anchor', 'middle')
    .attr('font-size', '13px')
    .text('False Positive Rate');

  g.append('text')
    .attr('transform', 'rotate(-90)')
    .attr('x', -height / 2)
    .attr('y', -45)
    .attr('text-anchor', 'middle')
    .attr('font-size', '13px')
    .text('True Positive Rate');

  // Title
  if (options.title) {
    g.append('text')
      .attr('x', width / 2)
      .attr('y', -10)
      .attr('text-anchor', 'middle')
      .attr('font-size', '14px')
      .attr('font-weight', 'bold')
      .text(options.title);
  }

  // Diagonal reference line (random classifier)
  g.append('line')
    .attr('x1', xScale(0))
    .attr('y1', yScale(0))
    .attr('x2', xScale(1))
    .attr('y2', yScale(1))
    .attr('stroke', '#999')
    .attr('stroke-dasharray', '6,4')
    .attr('stroke-width', 1);

  // Build point data
  const points = fpr.map((f, i) => ({
    fpr: f,
    tpr: tpr[i],
    threshold: thresholds[i] || 0,
  }));

  // Shaded AUC area
  const areaGen = d3.area()
    .x(d => xScale(d.fpr))
    .y0(height)
    .y1(d => yScale(d.tpr));

  g.append('path')
    .datum(points)
    .attr('fill', '#2196f3')
    .attr('fill-opacity', 0.15)
    .attr('d', areaGen);

  // ROC curve line
  const lineGen = d3.line()
    .x(d => xScale(d.fpr))
    .y(d => yScale(d.tpr));

  g.append('path')
    .datum(points)
    .attr('fill', 'none')
    .attr('stroke', '#1565c0')
    .attr('stroke-width', 2.5)
    .attr('d', lineGen);

  // AUC text
  if (auc !== undefined && auc !== null) {
    g.append('text')
      .attr('x', width - 10)
      .attr('y', height - 20)
      .attr('text-anchor', 'end')
      .attr('font-size', '14px')
      .attr('font-weight', 'bold')
      .attr('fill', '#1565c0')
      .text(`AUC = ${auc.toFixed(4)}`);
  }

  // Recommended threshold marker
  if (options.recommendedThreshold !== undefined) {
    const recThresh = options.recommendedThreshold;
    // Find closest point
    let bestIdx = 0;
    let bestDist = Infinity;
    for (let i = 0; i < thresholds.length; i++) {
      const dist = Math.abs(thresholds[i] - recThresh);
      if (dist < bestDist) {
        bestDist = dist;
        bestIdx = i;
      }
    }
    const recFPR = fpr[bestIdx];
    const recTPR = tpr[bestIdx];

    // Draw marker
    g.append('circle')
      .attr('cx', xScale(recFPR))
      .attr('cy', yScale(recTPR))
      .attr('r', 8)
      .attr('fill', '#e53935')
      .attr('fill-opacity', 0.8)
      .attr('stroke', '#fff')
      .attr('stroke-width', 2);

    g.append('text')
      .attr('x', xScale(recFPR) + 12)
      .attr('y', yScale(recTPR) + 4)
      .attr('font-size', '11px')
      .attr('fill', '#e53935')
      .text(`t=${recThresh.toFixed(3)}`);
  }

  // Tooltip for clickable points
  const tooltip = d3.select(container)
    .append('div')
    .style('position', 'absolute')
    .style('background', 'rgba(0,0,0,0.8)')
    .style('color', '#fff')
    .style('padding', '6px 10px')
    .style('border-radius', '4px')
    .style('font-size', '12px')
    .style('pointer-events', 'none')
    .style('opacity', 0);

  // Clickable scatter points (every 5th to avoid overplotting)
  const step = Math.max(1, Math.floor(points.length / 30));
  const sparsePoints = points.filter((_, i) => i % step === 0);

  g.selectAll('.roc-point')
    .data(sparsePoints)
    .join('circle')
    .attr('class', 'roc-point')
    .attr('cx', d => xScale(d.fpr))
    .attr('cy', d => yScale(d.tpr))
    .attr('r', 4)
    .attr('fill', '#1565c0')
    .attr('fill-opacity', 0.6)
    .attr('stroke', '#fff')
    .attr('stroke-width', 1)
    .attr('cursor', 'pointer')
    .on('mouseover', function(event, d) {
      d3.select(this).attr('r', 6).attr('fill-opacity', 1);
      tooltip
        .style('opacity', 1)
        .html(`Threshold: ${d.threshold.toFixed(4)}<br>FPR: ${d.fpr.toFixed(4)}<br>TPR: ${d.tpr.toFixed(4)}`);
    })
    .on('mousemove', function(event) {
      tooltip
        .style('left', (event.offsetX + 15) + 'px')
        .style('top', (event.offsetY - 10) + 'px');
    })
    .on('mouseout', function() {
      d3.select(this).attr('r', 4).attr('fill-opacity', 0.6);
      tooltip.style('opacity', 0);
    })
    .on('click', function(event, d) {
      if (options.onPointClick) {
        options.onPointClick(d.threshold);
      }
    });
}


/**
 * Render a recommendation table for threshold optimization results.
 *
 * recommendations: array of { column, threshold, sensitivity, specificity,
 *                              f1, fpr, true_positives, false_positives,
 *                              reads_affected }
 */
function renderThresholdTable(containerId, recommendations) {
  const container = document.getElementById(containerId);
  if (!container) return;
  container.innerHTML = '';

  if (!recommendations || recommendations.length === 0) {
    container.innerHTML = '<p>No threshold recommendations available.</p>';
    return;
  }

  let html = `
    <table>
      <caption>Recommended Thresholds</caption>
      <thead>
        <tr>
          <th>Confidence Column</th>
          <th>Threshold</th>
          <th>Sensitivity</th>
          <th>Specificity</th>
          <th>F1 Score</th>
          <th>FPR</th>
          <th>True Pos.</th>
          <th>False Pos.</th>
          <th>Reads Affected</th>
        </tr>
      </thead>
      <tbody>
  `;

  for (const rec of recommendations) {
    const f1Class = rec.f1 >= 0.9 ? 'badge-ok' : '';
    html += `
      <tr>
        <td><code>${rec.column || 'N/A'}</code></td>
        <td><strong>${rec.threshold.toFixed(4)}</strong></td>
        <td>${rec.sensitivity.toFixed(4)}</td>
        <td>${rec.specificity.toFixed(4)}</td>
        <td class="${f1Class}">${rec.f1.toFixed(4)}</td>
        <td>${rec.fpr.toFixed(4)}</td>
        <td>${(rec.true_positives || 0).toLocaleString()}</td>
        <td>${(rec.false_positives || 0).toLocaleString()}</td>
        <td>${(rec.reads_affected || 0).toLocaleString()}</td>
      </tr>
    `;
  }

  html += '</tbody></table>';
  container.innerHTML = html;
}
