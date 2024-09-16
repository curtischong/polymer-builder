import { atomicNumberToHexColor } from "@/components/CrystalViewer/atomColors";
import { atomicRadii } from "@/components/CrystalViewer/atomSizes";
import {
    BufferGeometry,
    Line,
    LineBasicMaterial,
    Mesh,
    MeshStandardMaterial,
    Scene,
    SphereGeometry,
    Vector3,
} from "three";

// Function to create a material for a sphere with a given radius
function createSphereMaterial(radius: number, color: number): Mesh {
    const sphereGeometry = new SphereGeometry(radius, 32, 32);
    const sphereMaterial = new MeshStandardMaterial({
        color,
        metalness: 0.1, // Increase metalness for a shiny appearance
        roughness: 0.5, // Decrease roughness for a smoother surface
    });
    return new Mesh(sphereGeometry, sphereMaterial);
}
// Initialize a Map to store the sphere materials
export const spheresMap = new Map<number, Mesh>();

// Loop through the atomicRadii and create sphere materials
for (const [atomicNumberStr, radius] of Object.entries(atomicRadii)) {
    const atomicNumber = parseInt(atomicNumberStr);
    const color = atomicNumberToHexColor[atomicNumber];
    const sphereMaterial = createSphereMaterial(radius / 400, color);
    spheresMap.set(atomicNumber, sphereMaterial);
}

export const subtract = (vec1: number[], vec2: number[]): number[] => {
    return vec1.map((x, i) => x - vec2[i]);
};

export const norm = (vec: number[]): number => {
    return Math.sqrt(vec.reduce((acc, x) => acc + x * x, 0));
};

export const fracCoordToCartesianCoord = (
    fracCoord: number[],
    lattice: number[][],
): number[] => {
    const [lattice0, lattice1, lattice2] = lattice;

    return [
        fracCoord[0] * lattice0[0] +
            fracCoord[1] * lattice1[0] +
            fracCoord[2] * lattice2[0],
        fracCoord[0] * lattice0[1] +
            fracCoord[1] * lattice1[1] +
            fracCoord[2] * lattice2[1],
        fracCoord[0] * lattice0[2] +
            fracCoord[1] * lattice1[2] +
            fracCoord[2] * lattice2[2],
    ];
};

export const findMiddleOfCoords = (coords: number[][]): number[] => {
    const smallest: number[] = [
        Number.POSITIVE_INFINITY,
        Number.POSITIVE_INFINITY,
        Number.POSITIVE_INFINITY,
    ];
    const largest: number[] = [
        Number.NEGATIVE_INFINITY,
        Number.NEGATIVE_INFINITY,
        Number.NEGATIVE_INFINITY,
    ];
    for (let i = 0; i < coords.length; i++) {
        const coord = coords[i];
        for (let j = 0; j < 3; j++) {
            const val = coord[j];
            if (val < smallest[j]) {
                smallest[j] = val;
            }
            if (val > largest[j]) {
                largest[j] = val;
            }
        }
    }
    return [
        (smallest[0] + largest[0]) / 2,
        (smallest[1] + largest[1]) / 2,
        (smallest[2] + largest[2]) / 2,
    ];
};

export const getParallelepipedCornerCoords = (
    lattice: number[][],
): number[][] => {
    const v1: number[] = lattice[0];
    const v2: number[] = lattice[1];
    const v3: number[] = lattice[2];

    const cornerCoords = [
        [0, 0, 0],
        v1,
        v1.map((val, idx) => val + v2[idx]),
        v2,
        v3,
        v1.map((val, idx) => val + v3[idx]),
        v1.map((val, idx) => val + v2[idx] + v3[idx]),
        v2.map((val, idx) => val + v3[idx]),
    ];
    return cornerCoords;
};

// Type definition for an edge
export type Edge = [number[], number[]];
export const getParallelepipedLatticeEdges = (points: number[][]): Edge[] => {
    // Create the edges of the parallelepiped as tuples of Cartesian coordinates
    let edges: Edge[] = [
        [points[0], points[1]],
        [points[1], points[2]],
        [points[2], points[3]],
        [points[3], points[0]],
        [points[4], points[5]],
        [points[5], points[6]],
        [points[6], points[7]],
        [points[7], points[4]],
        [points[0], points[4]],
        [points[1], points[5]],
        [points[2], points[6]],
        [points[3], points[7]],
    ];
    return edges;
};

export const drawEdge = (scene: Scene, edge: Edge) => {
    const p1 = edge[0];
    const p2 = edge[1];
    // TODO: create the vector 3 earlier so we don't need to create it here when drawing the line
    const points = [
        new Vector3(p1[0], p1[1], p1[2]),
        new Vector3(p2[0], p2[1], p2[2]),
    ];
    const geometry = new BufferGeometry().setFromPoints(points);

    // Create the material for the line
    const material = new LineBasicMaterial({ color: 0x0000ff });
    // Create the line and add it to the scene
    const line = new Line(geometry, material);
    scene.add(line);
};
