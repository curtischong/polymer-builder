import React from "react";
import { LineChart } from "@mantine/charts";
import { Paper, Text } from "@mantine/core";
import "@mantine/charts/styles.css";

type Force = number[]; // [x, y, z]
type ForceSet = Force[];

interface ForceLineChartProps {
    forceSets: ForceSet[];
}

interface ChartDataPoint {
    index: number;
    sum: number;
}

const ForceLineChart: React.FC<ForceLineChartProps> = ({ forceSets }) => {
    const data: ChartDataPoint[] = forceSets.map((forceSet, index) => ({
        index: index + 1,
        sum: forceSet.map((f) => Math.abs(f[2])).reduce((a, b) => a + b, 0),
    }));

    return (
        <Paper className="p-6 w-full max-w-3xl mx-auto rounded-lg shadow-md">
            <Text className="text-2xl font-bold mb-6 text-center text-gray-800">
                Sum of Absolute Z-Index Forces
            </Text>
            <LineChart
                h={350}
                data={data}
                dataKey="index"
                series={[
                    {
                        name: "sum",
                        color: "blue",
                    },
                ]}
                curveType="natural"
                // withLegend
                // withTooltip
                withDots={true}
                yAxisProps={{
                    tickSize: 1,
                    axisLine: { stroke: "#888888" },
                    ticks: [],
                }}
                xAxisProps={{
                    tickSize: 1,
                    axisLine: { stroke: "#888888" },
                    ticks: [],
                }}
                gridAxis="none"
                tooltipProps={{
                    content: ({ payload }) => {
                        if (payload && payload.length > 0) {
                            const dataPoint = payload[0]
                                .payload as ChartDataPoint;
                            return (
                                <div
                                    style={{
                                        background: "white",
                                        padding: "8px",
                                        borderRadius: "4px",
                                        boxShadow: "0 2px 4px rgba(0,0,0,0.1)",
                                        color: "#444444",
                                    }}
                                >
                                    <p style={{ fontWeight: "bold" }}>
                                        Index: {dataPoint.index}
                                    </p>
                                    <p>Sum: {dataPoint.sum.toFixed(2)}</p>
                                </div>
                            );
                        }
                        return null;
                    },
                }}
            />
        </Paper>
    );
};

export default ForceLineChart;
