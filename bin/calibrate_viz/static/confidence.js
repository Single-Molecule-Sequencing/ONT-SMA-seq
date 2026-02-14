/**
 * D3.js confidence visualization with draggable threshold lines.
 *
 * Provides:
 *   renderConfidenceScatter(containerId, data, options) - scatter of bc_start_conf vs bc_end_conf
 *   renderConfidenceKDE(containerId, data, options) - KDE of confidence scores
 *
 * Data format (scatter):
 *   { reads: [{ bc_start_conf, bc_end_conf, trunc_level }, ...] }
 *
 * Options (scatter):
 *   { xlabel, ylabel, startThreshold, flThreshold,
 *     onThresholdChange(startMin, flThresh) }
 *
 * Data format (KDE):
 *   { groups: { "label": { x: [...], y: [...], count: N }, ... },
 *     peaks: { "label": [...], ... } }
 *
 * Options (KDE):
 *   { xlabel, ylabel, thresholds: [{ value, label, color }] }
 */

// Color palette for trunc_level categories
const TRUNC_COLORS = {
  full_length: '#4caf50',
  bc1_target_bc2: '#2196f3',
  bc1_target: '#ff9800',
  adapter_only: '#e53935',
  unknown: '#9e9e9e',
};


/**
 * Render a scatter plot of bc_start_conf vs bc_end_conf colored by trunc_level.
 * Includes draggable horizontal and vertical threshold lines.
 */
function renderConfidenceScatter(containerId, data, options) {
  const container = document.getElementById(containerId);
  if (!container) return;
  container.innerHTML = '';

  const reads = data.reads || [];
  if (reads.length === 0) {
    container.innerHTML = '<p>No confidence data available.</p>';
    return;
  }

  const margin = { top: 20, right: 160, bottom: 50, left: 60 };
  const width = Math.max(container.clientWidth || 700, 500) - margin.left - margin.right;
  const height = 400 - margin.top - margin.bottom;

  const xScale = d3.scaleLinear().domain([0, 1]).range([0, width]);
  const yScale = d3.scaleLinear().domain([0, 1]).range([height, 0]);

  const svg = d3.select(container)
    .append('svg')
    .attr('width', width + margin.left + margin.right)
    .attr('height', height + margin.top + margin.bottom);

  const g = svg.append('g')
    .attr('transform', `translate(${margin.left},${margin.top})`);

  // Clip path to keep dots inside axes
  g.append('defs').append('clipPath')
    .attr('id', 'scatter-clip')
    .append('rect')
    .attr('width', width)
    .attr('height', height);

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
    .text(options.xlabel || 'bc_start_conf');

  g.append('text')
    .attr('transform', 'rotate(-90)')
    .attr('x', -height / 2)
    .attr('y', -45)
    .attr('text-anchor', 'middle')
    .attr('font-size', '13px')
    .text(options.ylabel || 'bc_end_conf');

  // Draw scatter points
  const plotArea = g.append('g').attr('clip-path', 'url(#scatter-clip)');

  plotArea.selectAll('.conf-dot')
    .data(reads)
    .join('circle')
    .attr('class', 'conf-dot')
    .attr('cx', d => xScale(d.bc_start_conf))
    .attr('cy', d => yScale(d.bc_end_conf != null ? d.bc_end_conf : 0))
    .attr('r', 3)
    .attr('fill', d => TRUNC_COLORS[d.trunc_level] || TRUNC_COLORS.unknown)
    .attr('opacity', 0.6);

  // Threshold state
  let startThresh = options.startThreshold || 0.6;
  let flThresh = options.flThreshold || 0.75;

  // --- Vertical threshold line (start_barcode_min) ---
  const vGroup = g.append('g').attr('class', 'threshold-group');

  const vLine = vGroup.append('line')
    .attr('class', 'threshold-line')
    .attr('x1', xScale(startThresh))
    .attr('x2', xScale(startThresh))
    .attr('y1', 0)
    .attr('y2', height)
    .attr('stroke', '#e53935');

  const vLabel = vGroup.append('text')
    .attr('x', xScale(startThresh) + 5)
    .attr('y', height - 5)
    .attr('font-size', '11px')
    .attr('fill', '#e53935')
    .text(`start_min: ${startThresh.toFixed(2)}`);

  const vHandle = vGroup.append('polygon')
    .attr('class', 'threshold-handle')
    .attr('points', _trianglePoints(xScale(startThresh), 0, 8))
    .attr('fill', '#e53935')
    .attr('opacity', 0.7);

  const vDragRect = vGroup.append('rect')
    .attr('x', xScale(startThresh) - 6)
    .attr('y', 0)
    .attr('width', 12)
    .attr('height', height)
    .attr('fill', 'transparent')
    .attr('cursor', 'ew-resize');

  vDragRect.call(d3.drag()
    .on('drag', function(event) {
      const newX = Math.max(0, Math.min(width, event.x));
      const newVal = xScale.invert(newX);
      vLine.attr('x1', newX).attr('x2', newX);
      vDragRect.attr('x', newX - 6);
      vHandle.attr('points', _trianglePoints(newX, 0, 8));
      vLabel.attr('x', newX + 5).text(`start_min: ${newVal.toFixed(2)}`);
    })
    .on('end', function(event) {
      const newX = Math.max(0, Math.min(width, event.x));
      startThresh = xScale.invert(newX);
      if (options.onThresholdChange) {
        options.onThresholdChange(startThresh, flThresh);
      }
    })
  );

  // --- Horizontal threshold line (full_length_threshold) ---
  const hGroup = g.append('g').attr('class', 'threshold-group');

  const hLine = hGroup.append('line')
    .attr('class', 'threshold-line')
    .attr('x1', 0)
    .attr('x2', width)
    .attr('y1', yScale(flThresh))
    .attr('y2', yScale(flThresh))
    .attr('stroke', '#1565c0');

  const hLabel = hGroup.append('text')
    .attr('x', 5)
    .attr('y', yScale(flThresh) - 5)
    .attr('font-size', '11px')
    .attr('fill', '#1565c0')
    .text(`fl_thresh: ${flThresh.toFixed(2)}`);

  const hHandle = hGroup.append('polygon')
    .attr('class', 'threshold-handle')
    .attr('points', _hTrianglePoints(0, yScale(flThresh), 8))
    .attr('fill', '#1565c0')
    .attr('opacity', 0.7);

  const hDragRect = hGroup.append('rect')
    .attr('x', 0)
    .attr('y', yScale(flThresh) - 6)
    .attr('width', width)
    .attr('height', 12)
    .attr('fill', 'transparent')
    .attr('cursor', 'ns-resize');

  hDragRect.call(d3.drag()
    .on('drag', function(event) {
      const newY = Math.max(0, Math.min(height, event.y));
      const newVal = yScale.invert(newY);
      hLine.attr('y1', newY).attr('y2', newY);
      hDragRect.attr('y', newY - 6);
      hHandle.attr('points', _hTrianglePoints(0, newY, 8));
      hLabel.attr('y', newY - 5).text(`fl_thresh: ${newVal.toFixed(2)}`);
    })
    .on('end', function(event) {
      const newY = Math.max(0, Math.min(height, event.y));
      flThresh = yScale.invert(newY);
      if (options.onThresholdChange) {
        options.onThresholdChange(startThresh, flThresh);
      }
    })
  );

  // Legend
  const legend = g.append('g')
    .attr('class', 'kde-legend')
    .attr('transform', `translate(${width + 10}, 10)`);

  const levels = ['full_length', 'bc1_target_bc2', 'bc1_target', 'adapter_only'];
  levels.forEach((level, i) => {
    const lg = legend.append('g')
      .attr('transform', `translate(0, ${i * 20})`);

    lg.append('rect')
      .attr('width', 14)
      .attr('height', 14)
      .attr('fill', TRUNC_COLORS[level]);

    lg.append('text')
      .attr('x', 20)
      .attr('y', 11)
      .attr('font-size', '12px')
      .text(level);
  });
}


/**
 * Render a KDE plot for confidence scores (bc_start_conf or bc_end_conf).
 * Supports multiple groups colored by trunc_level or assignment correctness.
 */
function renderConfidenceKDE(containerId, data, options) {
  const container = document.getElementById(containerId);
  if (!container) return;
  container.innerHTML = '';

  const groups = data.groups || {};
  const peaks = data.peaks || {};
  const groupNames = Object.keys(groups);
  if (groupNames.length === 0) {
    container.innerHTML = '<p>No data available.</p>';
    return;
  }

  const margin = { top: 20, right: 160, bottom: 50, left: 60 };
  const width = Math.max(container.clientWidth || 700, 500) - margin.left - margin.right;
  const height = 250 - margin.top - margin.bottom;

  // Use trunc_level colors if group names match, otherwise Tableau10
  const colorFn = (name) => TRUNC_COLORS[name] || d3.schemeTableau10[groupNames.indexOf(name) % 10];

  // Compute global x/y extent
  let xMin = Infinity, xMax = -Infinity, yMax = 0;
  for (const name of groupNames) {
    const g = groups[name];
    if (g.x.length === 0) continue;
    xMin = Math.min(xMin, d3.min(g.x));
    xMax = Math.max(xMax, d3.max(g.x));
    yMax = Math.max(yMax, d3.max(g.y));
  }
  if (!isFinite(xMin)) return;

  const xScale = d3.scaleLinear().domain([xMin, xMax]).range([0, width]);
  const yScale = d3.scaleLinear().domain([0, yMax * 1.1]).range([height, 0]);

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
    .call(d3.axisBottom(xScale).ticks(8));

  g.append('g')
    .attr('class', 'kde-axis')
    .call(d3.axisLeft(yScale).ticks(5));

  // Axis labels
  g.append('text')
    .attr('x', width / 2)
    .attr('y', height + 40)
    .attr('text-anchor', 'middle')
    .attr('font-size', '13px')
    .text(options.xlabel || 'Confidence');

  g.append('text')
    .attr('transform', 'rotate(-90)')
    .attr('x', -height / 2)
    .attr('y', -45)
    .attr('text-anchor', 'middle')
    .attr('font-size', '13px')
    .text(options.ylabel || 'Density');

  // Draw KDE curves
  for (const name of groupNames) {
    const gData = groups[name];
    const points = gData.x.map((xv, i) => ({ x: xv, y: gData.y[i] }));

    g.append('path')
      .datum(points)
      .attr('class', 'kde-curve')
      .attr('stroke', colorFn(name))
      .attr('d', d3.line()
        .x(d => xScale(d.x))
        .y(d => yScale(d.y))
      );

    // Peak markers
    if (peaks[name]) {
      for (const peakX of peaks[name]) {
        if (peakX >= xMin && peakX <= xMax) {
          g.append('line')
            .attr('class', 'peak-marker')
            .attr('x1', xScale(peakX))
            .attr('x2', xScale(peakX))
            .attr('y1', 0)
            .attr('y2', height)
            .attr('stroke', colorFn(name))
            .attr('opacity', 0.4);
        }
      }
    }
  }

  // Threshold lines (non-draggable, just markers in KDE)
  if (options.thresholds) {
    for (const th of options.thresholds) {
      if (th.value >= xMin && th.value <= xMax) {
        g.append('line')
          .attr('class', 'threshold-line')
          .attr('x1', xScale(th.value))
          .attr('x2', xScale(th.value))
          .attr('y1', 0)
          .attr('y2', height)
          .attr('stroke', th.color || '#e53935')
          .attr('stroke-dasharray', '4,4');

        g.append('text')
          .attr('x', xScale(th.value) + 4)
          .attr('y', 12)
          .attr('font-size', '11px')
          .attr('fill', th.color || '#e53935')
          .text(th.label || th.value.toFixed(2));
      }
    }
  }

  // Legend
  const legend = g.append('g')
    .attr('class', 'kde-legend')
    .attr('transform', `translate(${width + 10}, 10)`);

  groupNames.forEach((name, i) => {
    const lg = legend.append('g')
      .attr('transform', `translate(0, ${i * 20})`);

    lg.append('rect')
      .attr('width', 14)
      .attr('height', 14)
      .attr('fill', colorFn(name));

    const count = groups[name].count || 0;
    lg.append('text')
      .attr('x', 20)
      .attr('y', 11)
      .attr('font-size', '12px')
      .text(`${name} (${count.toLocaleString()})`);
  });
}


/**
 * Helper: generate triangle points for a vertical drag handle (pointing up).
 */
function _hTrianglePoints(cx, cy, size) {
  return `${cx},${cy} ${cx - size * 1.5},${cy - size} ${cx - size * 1.5},${cy + size}`;
}
