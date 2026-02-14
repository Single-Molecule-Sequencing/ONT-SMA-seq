/**
 * D3.js KDE distribution visualization with draggable thresholds.
 *
 * Data format:
 *   { groups: { "nb05": { x: [...], y: [...], count: N }, ... },
 *     peaks:  { "nb05": [123, 456], ... } }
 *
 * Options:
 *   { xlabel, ylabel, column, groupBy,
 *     thresholds: [{ value, label, color, statsTarget }],
 *     expectedSizes: [{ value, label }] }
 */

function renderKDE(containerId, data, options) {
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

  // Dimensions
  const margin = { top: 20, right: 180, bottom: 50, left: 60 };
  const width = Math.max(container.clientWidth || 800, 600) - margin.left - margin.right;
  const height = 320 - margin.top - margin.bottom;

  const color = d3.scaleOrdinal(d3.schemeTableau10).domain(groupNames);

  // Compute global x/y extent
  let xMin = Infinity, xMax = -Infinity, yMax = 0;
  for (const name of groupNames) {
    const g = groups[name];
    const xs = g.x, ys = g.y;
    if (xs.length === 0) continue;
    xMin = Math.min(xMin, d3.min(xs));
    xMax = Math.max(xMax, d3.max(xs));
    yMax = Math.max(yMax, d3.max(ys));
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
    .text(options.xlabel || '');

  g.append('text')
    .attr('transform', 'rotate(-90)')
    .attr('x', -height / 2)
    .attr('y', -45)
    .attr('text-anchor', 'middle')
    .attr('font-size', '13px')
    .text(options.ylabel || 'Density');

  // Line generator
  const line = d3.line()
    .x((d, i) => xScale(groups[d.__group].x[i]))
    .y((d, i) => yScale(groups[d.__group].y[i]));

  // Draw KDE curves
  for (const name of groupNames) {
    const gData = groups[name];
    const points = gData.x.map((xv, i) => ({ x: xv, y: gData.y[i] }));

    g.append('path')
      .datum(points)
      .attr('class', 'kde-curve')
      .attr('stroke', color(name))
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
            .attr('stroke', color(name))
            .attr('opacity', 0.4);
        }
      }
    }
  }

  // Expected size markers
  if (options.expectedSizes) {
    for (const es of options.expectedSizes) {
      if (es.value >= xMin && es.value <= xMax) {
        g.append('line')
          .attr('class', 'expected-marker')
          .attr('x1', xScale(es.value))
          .attr('x2', xScale(es.value))
          .attr('y1', 0)
          .attr('y2', height);

        g.append('text')
          .attr('x', xScale(es.value) + 4)
          .attr('y', 12)
          .attr('font-size', '11px')
          .attr('fill', '#999')
          .text(es.label || es.value);
      }
    }
  }

  // Draggable threshold lines
  if (options.thresholds) {
    for (const th of options.thresholds) {
      const thGroup = g.append('g')
        .attr('class', 'threshold-group');

      const thLine = thGroup.append('line')
        .attr('class', 'threshold-line')
        .attr('x1', xScale(th.value))
        .attr('x2', xScale(th.value))
        .attr('y1', 0)
        .attr('y2', height)
        .attr('stroke', th.color || '#e53935');

      // Drag handle (triangle/diamond at top)
      const handleSize = 8;
      thGroup.append('polygon')
        .attr('class', 'threshold-handle')
        .attr('points', _trianglePoints(xScale(th.value), 0, handleSize))
        .attr('stroke', th.color || '#e53935')
        .attr('fill', th.color || '#e53935')
        .attr('opacity', 0.7);

      // Label
      const thLabel = thGroup.append('text')
        .attr('x', xScale(th.value) + 5)
        .attr('y', height - 5)
        .attr('font-size', '11px')
        .attr('fill', th.color || '#e53935')
        .text(th.label ? `${th.label}: ${th.value.toFixed(0)}` : th.value.toFixed(0));

      // Invisible wide rect for easier dragging
      const dragRect = thGroup.append('rect')
        .attr('x', xScale(th.value) - 6)
        .attr('y', 0)
        .attr('width', 12)
        .attr('height', height)
        .attr('fill', 'transparent')
        .attr('cursor', 'ew-resize');

      // Drag behavior
      const drag = d3.drag()
        .on('drag', function(event) {
          let newX = Math.max(0, Math.min(width, event.x));
          let newVal = xScale.invert(newX);

          thLine.attr('x1', newX).attr('x2', newX);
          dragRect.attr('x', newX - 6);
          thGroup.select('polygon')
            .attr('points', _trianglePoints(newX, 0, handleSize));
          thLabel
            .attr('x', newX + 5)
            .text(th.label ? `${th.label}: ${newVal.toFixed(0)}` : newVal.toFixed(0));
        })
        .on('end', function(event) {
          let newX = Math.max(0, Math.min(width, event.x));
          let newVal = xScale.invert(newX);
          th.value = newVal;

          // Fire HTMX request for threshold impact
          if (th.statsTarget) {
            const params = new URLSearchParams({
              column: options.column || 'readlen',
              threshold: newVal.toFixed(2),
              db: window.currentDb || 0,
            });
            if (options.groupBy) {
              params.set('group_by', options.groupBy);
            }

            fetch('/api/threshold-impact?' + params)
              .then(r => r.text())
              .then(html => {
                const target = document.getElementById(th.statsTarget);
                if (target) target.innerHTML = html;
              });
          }
        });

      dragRect.call(drag);
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
      .attr('fill', color(name));

    const count = groups[name].count || 0;
    lg.append('text')
      .attr('x', 20)
      .attr('y', 11)
      .attr('font-size', '12px')
      .text(`${name} (${count.toLocaleString()})`);
  });
}


/**
 * Render end-reason stacked area subplot sharing an x-axis reference.
 */
function renderEndReasonSubplot(containerId, data, xScaleRef) {
  const container = document.getElementById(containerId);
  if (!container) return;
  container.innerHTML = '';

  const groups = data.groups || {};
  const groupNames = Object.keys(groups);
  if (groupNames.length === 0) return;

  const margin = { top: 10, right: 180, bottom: 30, left: 60 };
  const width = Math.max(container.clientWidth || 800, 600) - margin.left - margin.right;
  const height = 140 - margin.top - margin.bottom;

  // End-reason color palette
  const erColors = {
    'signal_positive': '#4caf50',
    'data_service_unblock_mux_change': '#f44336',
    'unblock_mux_change': '#ff9800',
    'signal_negative': '#9e9e9e',
    'mux_change': '#2196f3',
  };
  const color = d3.scaleOrdinal()
    .domain(groupNames)
    .range(groupNames.map(n => erColors[n] || d3.schemeTableau10[groupNames.indexOf(n) % 10]));

  // Compute global x extent
  let xMin = Infinity, xMax = -Infinity;
  for (const name of groupNames) {
    const xs = groups[name].x;
    if (xs.length === 0) continue;
    xMin = Math.min(xMin, d3.min(xs));
    xMax = Math.max(xMax, d3.max(xs));
  }
  if (!isFinite(xMin)) return;

  const xScale = xScaleRef || d3.scaleLinear().domain([xMin, xMax]).range([0, width]);

  // Resample all groups to a common x grid for stacking
  const nPoints = 256;
  const xGrid = d3.range(nPoints).map(i => xMin + (xMax - xMin) * i / (nPoints - 1));

  // Interpolate each group onto the grid
  const interpolated = {};
  for (const name of groupNames) {
    const g = groups[name];
    const interp = [];
    for (const xv of xGrid) {
      // Find nearest pair in g.x and linearly interpolate g.y
      let val = 0;
      const xs = g.x, ys = g.y;
      if (xs.length === 0) { interp.push(0); continue; }
      if (xv <= xs[0]) { val = ys[0]; }
      else if (xv >= xs[xs.length - 1]) { val = ys[ys.length - 1]; }
      else {
        for (let j = 0; j < xs.length - 1; j++) {
          if (xv >= xs[j] && xv <= xs[j + 1]) {
            const t = (xv - xs[j]) / (xs[j + 1] - xs[j]);
            val = ys[j] + t * (ys[j + 1] - ys[j]);
            break;
          }
        }
      }
      // Weight by count for stacking
      interp.push(val * (g.count || 1));
    }
    interpolated[name] = interp;
  }

  // Build stack data
  const stackData = xGrid.map((xv, i) => {
    const row = { x: xv };
    for (const name of groupNames) {
      row[name] = interpolated[name][i];
    }
    return row;
  });

  const stack = d3.stack().keys(groupNames);
  const series = stack(stackData);

  // Y scale
  const yMax = d3.max(series, s => d3.max(s, d => d[1])) || 1;
  const yScale = d3.scaleLinear().domain([0, yMax]).range([height, 0]);

  const svg = d3.select(container)
    .append('svg')
    .attr('width', width + margin.left + margin.right)
    .attr('height', height + margin.top + margin.bottom);

  const gEl = svg.append('g')
    .attr('transform', `translate(${margin.left},${margin.top})`);

  // Axes
  gEl.append('g')
    .attr('class', 'kde-axis')
    .attr('transform', `translate(0,${height})`)
    .call(d3.axisBottom(xScale).ticks(8));

  gEl.append('g')
    .attr('class', 'kde-axis')
    .call(d3.axisLeft(yScale).ticks(3));

  // Label
  gEl.append('text')
    .attr('transform', 'rotate(-90)')
    .attr('x', -height / 2)
    .attr('y', -45)
    .attr('text-anchor', 'middle')
    .attr('font-size', '11px')
    .text('End Reason');

  // Stacked area
  const area = d3.area()
    .x((d, i) => xScale(stackData[i].x))
    .y0(d => yScale(d[0]))
    .y1(d => yScale(d[1]));

  gEl.selectAll('.er-area')
    .data(series)
    .join('path')
    .attr('class', 'er-area')
    .attr('fill', d => color(d.key))
    .attr('opacity', 0.7)
    .attr('d', area);

  // Legend
  const legend = gEl.append('g')
    .attr('class', 'kde-legend')
    .attr('transform', `translate(${width + 10}, 0)`);

  groupNames.forEach((name, i) => {
    const lg = legend.append('g')
      .attr('transform', `translate(0, ${i * 16})`);

    lg.append('rect')
      .attr('width', 12)
      .attr('height', 12)
      .attr('fill', color(name));

    lg.append('text')
      .attr('x', 16)
      .attr('y', 10)
      .attr('font-size', '10px')
      .text(name);
  });
}


/**
 * Helper: generate triangle points string for a drag handle.
 */
function _trianglePoints(cx, cy, size) {
  return `${cx},${cy} ${cx - size},${cy - size * 1.5} ${cx + size},${cy - size * 1.5}`;
}
