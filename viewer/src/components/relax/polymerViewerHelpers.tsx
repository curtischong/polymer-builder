export const getCornerCoords = (coords: number[][]): number[][] => {
    let minX = Number.POSITIVE_INFINITY;
    let maxX = Number.NEGATIVE_INFINITY;
    let minY = Number.POSITIVE_INFINITY;
    let maxY = Number.NEGATIVE_INFINITY;
    let minZ = Number.POSITIVE_INFINITY;
    let maxZ = Number.NEGATIVE_INFINITY;

    for (let i = 0; i < coords.length; i++) {
        const [x, y, z] = coords[i];

        minX = Math.min(minX, x);
        maxX = Math.max(maxX, x);

        minY = Math.min(minY, y);
        maxY = Math.max(maxY, y);

        minZ = Math.min(minZ, z);
        maxZ = Math.max(maxZ, z);
    }
    return [
        [minX, minY, minZ],
        [maxX, maxY, maxZ],
    ];
};
