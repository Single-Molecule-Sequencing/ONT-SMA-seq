/**
 * D3.js barcode separation visualization: pairwise edit distance heatmap
 * and sortable separation metrics table.
 *
 * Heatmap data format:
 *   { ids: ["nb01", "nb05", ...],
 *     matrix: { "nb01": { "nb01": 0, "nb05": 12, ... }, ... } }
 *
 * Metrics data format:
 *   { "nb01": { mean_true_ed, mean_next_best_ed, separation_gap,
 *               estimated_error_rate, n_reads }, ... }
 */

/**
 * Render a pairwise edit distance heatmap.
 *
 * @param {string} containerId - DOM element ID to render into.
 * @param {object} data - { ids: string[], matrix: {[id]: {[id]: number}} }
 */
function renderHeatmap(containerId, data) {
  const container = document.getElementById(containerId);
  if (!container) return;
  container.innerHTML = '';

  const ids = data.ids || [];
  const matrix = data.matrix || {};
  if (ids.length === 0) {
    container.innerHTML = '<p>No barcode data available.</p>';
    return;
  }

  // Find min/max distances (excluding diagonal zeros)
  let minDist = Infinity, maxDist = -Infinity;
  for (const a of ids) {
    for (const b of ids) {
      if (a === b) continue;
      const d = matrix[a][b];
      if (d < minDist) minDist = d;
      if (d > maxDist) maxDist = d;
    }
  }
  if (!isFinite(minDist)) { minDist = 0; maxDist = 1; }

  // Dimensions
  const cellSize = Math.max(12, Math.min(28, 600 / ids.length));
  const labelWidth = 50;
  const margin = { top: labelWidth + 10, right: 60, bottom: 20, left: labelWidth + 10 };
  const gridSize = cellSize * ids.length;
  const width = gridSize + margin.left + margin.right;
  const height = gridSize + margin.top + margin.bottom;

  // Color scale: low distance = red (danger), high distance = green (good separation)
  const colorScale = d3.scaleSequential()
    .domain([minDist, maxDist])
    .interpolator(d3.interpolateRgb('#e53935', '#43a047'));

  const svg = d3.select(container)
    .append('svg')
    .attr('width', width)
    .attr('height', height);

  const g = svg.append('g')
    .attr('transform', `translate(${margin.left},${margin.top})`);

  // Row labels (left)
  g.selectAll('.row-label')
    .data(ids)
    .join('text')
    .attr('class', 'row-label')
    .attr('x', -6)
    .attr('y', (d, i) => i * cellSize + cellSize / 2)
    .attr('text-anchor', 'end')
    .attr('dominant-baseline', 'middle')
    .attr('font-size', Math.min(11, cellSize - 2) + 'px')
    .text(d => d);

  // Column labels (top, rotated)
  g.selectAll('.col-label')
    .data(ids)
    .join('text')
    .attr('class', 'col-label')
    .attr('transform', (d, i) => `translate(${i * cellSize + cellSize / 2},-6) rotate(-45)`)
    .attr('text-anchor', 'start')
    .attr('font-size', Math.min(11, cellSize - 2) + 'px')
    .text(d => d);

  // Tooltip
  const tooltip = d3.select(container)
    .append('div')
    .style('position', 'absolute')
    .style('background', 'rgba(0,0,0,0.8)')
    .style('color', 'white')
    .style('padding', '4px 8px')
    .style('border-radius', '4px')
    .style('font-size', '12px')
    .style('pointer-events', 'none')
    .style('display', 'none')
    .style('z-index', '10');

  // Cells
  for (let i = 0; i < ids.length; i++) {
    for (let j = 0; j < ids.length; j++) {
      const a = ids[i];
      const b = ids[j];
      const dist = matrix[a][b];
      const isDiag = (i === j);

      g.append('rect')
        .attr('x', j * cellSize)
        .attr('y', i * cellSize)
        .attr('width', cellSize - 1)
        .attr('height', cellSize - 1)
        .attr('fill', isDiag ? '#e0e0e0' : colorScale(dist))
        .attr('stroke', '#fff')
        .attr('stroke-width', 0.5)
        .attr('rx', 1)
        .style('cursor', 'pointer')
        .on('mouseover', function(event) {
          tooltip
            .style('display', 'block')
            .html(`<strong>${a} vs ${b}</strong><br>Edit distance: ${dist}`);
        })
        .on('mousemove', function(event) {
          const rect = container.getBoundingClientRect();
          tooltip
            .style('left', (event.clientX - rect.left + 10) + 'px')
            .style('top', (event.clientY - rect.top - 28) + 'px');
        })
        .on('mouseout', function() {
          tooltip.style('display', 'none');
        })
        .on('click', function() {
          // Dispatch custom event for barcode detail panel
          container.dispatchEvent(new CustomEvent('barcode-click', {
            detail: { barcodeA: a, barcodeB: b, distance: dist },
            bubbles: true,
          }));
        });

      // Show distance text in cells if they're large enough
      if (cellSize >= 20 && !isDiag) {
        g.append('text')
          .attr('x', j * cellSize + cellSize / 2)
          .attr('y', i * cellSize + cellSize / 2)
          .attr('text-anchor', 'middle')
          .attr('dominant-baseline', 'middle')
          .attr('font-size', Math.min(9, cellSize - 6) + 'px')
          .attr('fill', dist < (minDist + maxDist) / 2 ? 'white' : '#333')
          .attr('pointer-events', 'none')
          .text(dist);
      }
    }
  }

  // Color legend
  const legendWidth = 15;
  const legendHeight = gridSize;
  const legendX = gridSize + 15;

  const legendScale = d3.scaleLinear()
    .domain([minDist, maxDist])
    .range([legendHeight, 0]);

  const legendAxis = d3.axisRight(legendScale).ticks(5).tickFormat(d3.format('d'));

  // Gradient
  const defs = svg.append('defs');
  const gradient = defs.append('linearGradient')
    .attr('id', 'heatmap-gradient')
    .attr('x1', '0%').attr('y1', '100%')
    .attr('x2', '0%').attr('y2', '0%');

  gradient.append('stop').attr('offset', '0%').attr('stop-color', '#e53935');
  gradient.append('stop').attr('offset', '100%').attr('stop-color', '#43a047');

  g.append('rect')
    .attr('x', legendX)
    .attr('y', 0)
    .attr('width', legendWidth)
    .attr('height', legendHeight)
    .style('fill', 'url(#heatmap-gradient)');

  g.append('g')
    .attr('transform', `translate(${legendX + legendWidth}, 0)`)
    .call(legendAxis)
    .selectAll('text')
    .attr('font-size', '10px');

  g.append('text')
    .attr('x', legendX + legendWidth / 2)
    .attr('y', -8)
    .attr('text-anchor', 'middle')
    .attr('font-size', '10px')
    .text('Edit dist.');

  // Make container position relative for tooltip
  d3.select(container).style('position', 'relative');
}


/**
 * Render a sortable separation metrics table.
 *
 * @param {string} containerId - DOM element ID to render into.
 * @param {object} metrics - { "nb01": { mean_true_ed, mean_next_best_ed,
 *                              separation_gap, estimated_error_rate, n_reads }, ... }
 */
function renderSeparationTable(containerId, metrics) {
  const container = document.getElementById(containerId);
  if (!container) return;
  container.innerHTML = '';

  const ids = Object.keys(metrics);
  if (ids.length === 0) {
    container.innerHTML = '<p>No separation metrics available.</p>';
    return;
  }

  // Build rows array for sorting
  let rows = ids.map(id => ({
    barcode: id,
    ...metrics[id],
  }));

  // State
  let sortCol = 'separation_gap';
  let sortAsc = true;

  function renderTable() {
    container.innerHTML = '';

    // Sort
    rows.sort((a, b) => {
      const va = a[sortCol], vb = b[sortCol];
      if (typeof va === 'string') return sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
      return sortAsc ? va - vb : vb - va;
    });

    const table = document.createElement('table');
    table.setAttribute('role', 'grid');

    // Header
    const thead = document.createElement('thead');
    const headerRow = document.createElement('tr');

    const columns = [
      { key: 'barcode', label: 'Barcode' },
      { key: 'n_reads', label: 'Reads' },
      { key: 'mean_true_ed', label: 'Mean Match ED' },
      { key: 'mean_next_best_ed', label: 'Mean Next-Best ED' },
      { key: 'separation_gap', label: 'Gap' },
      { key: 'estimated_error_rate', label: 'Est. Error Rate' },
    ];

    for (const col of columns) {
      const th = document.createElement('th');
      th.style.cursor = 'pointer';
      th.style.userSelect = 'none';

      let arrow = '';
      if (sortCol === col.key) {
        arrow = sortAsc ? ' \u25B2' : ' \u25BC';
      }
      th.textContent = col.label + arrow;

      th.addEventListener('click', () => {
        if (sortCol === col.key) {
          sortAsc = !sortAsc;
        } else {
          sortCol = col.key;
          sortAsc = true;
        }
        renderTable();
      });

      headerRow.appendChild(th);
    }
    thead.appendChild(headerRow);
    table.appendChild(thead);

    // Body
    const tbody = document.createElement('tbody');
    for (const row of rows) {
      const tr = document.createElement('tr');
      tr.style.cursor = 'pointer';

      // Highlight row based on error rate
      if (row.estimated_error_rate > 0.1) {
        tr.style.background = 'rgba(229, 57, 53, 0.08)';
      } else if (row.separation_gap < 5) {
        tr.style.background = 'rgba(255, 152, 0, 0.08)';
      }

      tr.addEventListener('click', () => {
        container.dispatchEvent(new CustomEvent('barcode-select', {
          detail: { barcode: row.barcode },
          bubbles: true,
        }));
      });

      const cells = [
        row.barcode,
        row.n_reads.toLocaleString(),
        row.mean_true_ed.toFixed(2),
        row.mean_next_best_ed.toFixed(2),
        row.separation_gap.toFixed(2),
        (row.estimated_error_rate * 100).toFixed(2) + '%',
      ];

      for (const val of cells) {
        const td = document.createElement('td');
        td.textContent = val;
        tr.appendChild(td);
      }

      tbody.appendChild(tr);
    }
    table.appendChild(tbody);
    container.appendChild(table);
  }

  renderTable();
}
