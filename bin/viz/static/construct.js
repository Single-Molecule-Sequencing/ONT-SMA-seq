/**
 * D3.js construct diagram for SMA-seq library product visualization.
 *
 * Provides three public functions:
 *   - renderConstruct(data, containerSel, detailSel)
 *   - renderTruncationLadder(data, containerSel)
 *   - showDetail(data, region, detailSel)
 */

/* global d3 */

// ---------------------------------------------------------------------------
// Region definitions
// ---------------------------------------------------------------------------

const REGIONS = [
  { key: "adapter_5",    label: "Adapter 5'",     color: "#607d8b", fraction: 0.10 },
  { key: "barcode1",     label: "Barcode 1",      color: "#2196f3", fraction: 0.10 },
  { key: "target",       label: "Target Insert",   color: "#4caf50", fraction: 0.50 },
  { key: "rc_barcode2",  label: "RC Barcode 2",   color: "#ff9800", fraction: 0.10 },
  { key: "rc_adapter_3", label: "RC Adapter 3'",  color: "#9e9e9e", fraction: 0.10 },
];

const GAP = 4;

// ---------------------------------------------------------------------------
// Truncation class definitions
// ---------------------------------------------------------------------------

const TRUNCATION_CLASSES = [
  { key: "full",          label: "Full",           weights: [1, 1, 1, 1, 1],       color: "#4caf50" },
  { key: "trunc_3prime",  label: "Trunc-3prime",   weights: [0, 1, 1, 1, 0],       color: "#8bc34a" },
  { key: "trunc_target",  label: "Trunc-target",   weights: [1, 1, 0.6, 0, 0],     color: "#ff9800" },
  { key: "trunc_barcode", label: "Trunc-barcode",  weights: [1, 1, 0, 0, 0],       color: "#ff5722" },
  { key: "adapter_only",  label: "Adapter-only",   weights: [1, 0, 0, 0, 0],       color: "#9e9e9e" },
  { key: "chimeric",      label: "Chimeric",       weights: [0.5, 0.5, 0.3, 0.5, 0], color: "#f44336" },
];

// ---------------------------------------------------------------------------
// renderConstruct
// ---------------------------------------------------------------------------

function renderConstruct(data, containerSel, detailSel) {
  const container = d3.select(containerSel);
  container.selectAll("*").remove();

  const width = 800;
  const height = 120;
  const margin = { top: 10, right: 10, bottom: 30, left: 10 };
  const drawWidth = width - margin.left - margin.right;
  const drawHeight = height - margin.top - margin.bottom;

  const svg = container.append("svg")
    .attr("viewBox", `0 0 ${width} ${height}`)
    .attr("preserveAspectRatio", "xMidYMid meet")
    .style("width", "100%")
    .style("max-width", `${width}px`);

  const g = svg.append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);

  // Compute total fraction for normalization
  const totalFraction = REGIONS.reduce((sum, r) => sum + r.fraction, 0);
  const totalGap = GAP * (REGIONS.length - 1);
  const usableWidth = drawWidth - totalGap;

  let x = 0;

  REGIONS.forEach((region, i) => {
    const w = (region.fraction / totalFraction) * usableWidth;
    const rectHeight = drawHeight * 0.65;
    const y = (drawHeight - rectHeight) / 2;

    // Rectangle
    g.append("rect")
      .attr("class", "construct-region")
      .attr("x", x)
      .attr("y", y)
      .attr("width", w)
      .attr("height", rectHeight)
      .attr("rx", 4)
      .attr("fill", region.color)
      .style("cursor", "pointer")
      .on("click", function () {
        showDetail(data, region.key, detailSel);
      });

    // Label text (white, centered inside rectangle)
    g.append("text")
      .attr("x", x + w / 2)
      .attr("y", y + rectHeight / 2)
      .attr("dy", "0.35em")
      .attr("text-anchor", "middle")
      .attr("fill", "white")
      .attr("font-size", w < 70 ? "9px" : "12px")
      .attr("font-weight", "600")
      .attr("pointer-events", "none")
      .text(region.label);

    x += w + GAP;
  });

  // Direction arrow below
  const arrowY = drawHeight + 5;
  const arrowG = g.append("g")
    .attr("transform", `translate(0,${arrowY})`);

  arrowG.append("text")
    .attr("x", 0)
    .attr("y", 0)
    .attr("dy", "0.8em")
    .attr("fill", "#666")
    .attr("font-size", "13px")
    .attr("font-family", "monospace")
    .text("5\u2032 \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2192 3\u2032");
}

// ---------------------------------------------------------------------------
// renderTruncationLadder
// ---------------------------------------------------------------------------

function renderTruncationLadder(data, containerSel) {
  const container = d3.select(containerSel);
  container.selectAll("*").remove();

  const width = 800;
  const rowHeight = 36;
  const labelWidth = 180;
  const margin = { top: 10, right: 10, bottom: 10, left: 10 };
  const barAreaWidth = width - margin.left - margin.right - labelWidth;
  const height = margin.top + margin.bottom + TRUNCATION_CLASSES.length * rowHeight;

  const svg = container.append("svg")
    .attr("viewBox", `0 0 ${width} ${height}`)
    .attr("preserveAspectRatio", "xMidYMid meet")
    .style("width", "100%")
    .style("max-width", `${width}px`);

  const g = svg.append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);

  // Compute total fraction for normalization (same as construct)
  const totalFraction = REGIONS.reduce((sum, r) => sum + r.fraction, 0);
  const totalGap = GAP * (REGIONS.length - 1);
  const usableBarWidth = barAreaWidth - totalGap;

  // Get truncation rules from data if available
  const rules = data && data.truncation && data.truncation.rules
    ? data.truncation.rules
    : {};

  TRUNCATION_CLASSES.forEach((tc, rowIdx) => {
    const y = rowIdx * rowHeight;
    const barHeight = rowHeight * 0.6;
    const barY = y + (rowHeight - barHeight) / 2;

    let x = 0;

    // Draw colored segments for each region
    REGIONS.forEach((region, i) => {
      const fullWidth = (region.fraction / totalFraction) * usableBarWidth;
      const weight = tc.weights[i];

      if (weight > 0) {
        const segWidth = fullWidth * Math.min(weight, 1);
        const opacity = weight < 1 ? 0.5 : 1.0;

        g.append("rect")
          .attr("x", x)
          .attr("y", barY)
          .attr("width", segWidth)
          .attr("height", barHeight)
          .attr("rx", 3)
          .attr("fill", region.color)
          .attr("opacity", opacity);
      }

      x += fullWidth + GAP;
    });

    // Class label on the right side
    const labelX = barAreaWidth + 8;
    g.append("text")
      .attr("x", labelX)
      .attr("y", barY + barHeight / 2)
      .attr("dy", "0.35em")
      .attr("fill", tc.color)
      .attr("font-size", "13px")
      .attr("font-weight", "600")
      .text(tc.label);

    // Assignment rule text below the label
    const ruleText = rules[tc.key] || "";
    if (ruleText) {
      g.append("text")
        .attr("x", labelX)
        .attr("y", barY + barHeight / 2 + 13)
        .attr("dy", "0.35em")
        .attr("fill", "#888")
        .attr("font-size", "10px")
        .text(ruleText);
    }
  });
}

// ---------------------------------------------------------------------------
// showDetail
// ---------------------------------------------------------------------------

function showDetail(data, region, detailSel) {
  const container = d3.select(detailSel);
  container.selectAll("*").remove();

  const construct = data && data.construct ? data.construct : {};
  const classification = data && data.classification ? data.classification : {};

  const article = container.append("article")
    .style("padding", "1rem")
    .style("border", "1px solid var(--pico-muted-border-color)")
    .style("border-radius", "0.5rem")
    .style("margin-top", "0.5rem");

  let title = "";
  let content = "";

  switch (region) {
    case "adapter_5":
      title = "Adapter 5\u2032";
      content = construct.adapter_5prime
        ? `<p><strong>Sequence:</strong> <code>${construct.adapter_5prime}</code></p>
           <p><strong>Length:</strong> ${construct.adapter_5prime.length} bp</p>`
        : "<p>No adapter 5\u2032 sequence configured.</p>";
      break;

    case "barcode1":
      title = "Barcode 1 (Start)";
      content = `<p><strong>Classification:</strong> Classified in first <strong>${classification.barcode_search_window || "N/A"}</strong> bp of the read.</p>
                 <p><strong>Confidence formula:</strong> <code>${classification.confidence_formula || "N/A"}</code></p>`;
      break;

    case "target":
      title = "Target Insert";
      content = `<p><strong>Insert type:</strong> ${construct.insert_type || "N/A"}</p>`;
      if (construct.notes) {
        content += `<p><strong>Notes:</strong> ${construct.notes}</p>`;
      }
      break;

    case "rc_barcode2":
      title = "RC Barcode 2 (End)";
      content = `<p><strong>Classification:</strong> Classified in last <strong>${classification.barcode_search_window || "N/A"}</strong> bp of the read.</p>
                 <p>This is the reverse-complement of the downstream barcode.</p>`;
      break;

    case "rc_adapter_3":
      title = "RC Adapter 3\u2032";
      content = construct.adapter_3prime
        ? `<p><strong>Sequence (original 3\u2032):</strong> <code>${construct.adapter_3prime}</code></p>
           <p><strong>Length:</strong> ${construct.adapter_3prime.length} bp</p>
           <p>Stored as reverse-complement in the read.</p>`
        : "<p>No adapter 3\u2032 sequence configured.</p>";
      break;

    default:
      title = "Unknown Region";
      content = "<p>No details available.</p>";
  }

  article.html(`<h4 style="margin-top:0">${title}</h4>${content}`);
}
