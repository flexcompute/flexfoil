/**
 * OutputPoints - Displays computed SSP output points along the s-axis
 */

import React from 'react';

interface OutputPointsProps {
  si: number[];
  toSvgX: (s: number) => number;
  yPosition: number;
}

const OUTPUT_POINT_RADIUS = 3;

export const OutputPoints: React.FC<OutputPointsProps> = ({ si, toSvgX, yPosition }) => {
  return (
    <g className="output-points">
      {/* Label */}
      <text
        x={toSvgX(0) - 35}
        y={yPosition + 4}
        fontSize="10"
        fill="var(--text-muted, #666)"
      >
        s<tspan fontSize="8" dy="2">i</tspan>
      </text>
      
      {/* Baseline */}
      <line
        x1={toSvgX(0)}
        y1={yPosition}
        x2={toSvgX(1)}
        y2={yPosition}
        stroke="var(--border-color, #ccc)"
        strokeWidth="1"
        strokeDasharray="2,2"
      />
      
      {/* Output points */}
      {si.map((s, i) => (
        <circle
          key={i}
          cx={toSvgX(s)}
          cy={yPosition}
          r={OUTPUT_POINT_RADIUS}
          fill="var(--accent-warning, #FF5722)"
          opacity={0.8}
        />
      ))}
    </g>
  );
};

export default OutputPoints;
