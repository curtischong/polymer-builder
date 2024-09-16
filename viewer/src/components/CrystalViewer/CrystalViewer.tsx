"use client";
import {
    drawEdge,
    Edge,
    findMiddleOfCoords,
    fracCoordToCartesianCoord,
    getParallelepipedCornerCoords,
    getParallelepipedLatticeEdges,
    norm,
    spheresMap,
    subtract,
} from "@/app/components/CrystalViewer/crystalViewerHelpers";
import React, { useRef, useEffect, useState } from "react";
import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";

export interface Frame {
    atomicNumbers: number[];
    fracCoord: number[][]; // PERF: maybe replace with a fixed size array if there's perf issues
    lattice: number[][];
}

export interface CrystalViewerProps {
    frames: Frame[];
}

// this function preprocesses the frames (e.g. visualize it as a mxmxm supercell / center all atoms),
// so when we pass them into the CrystalViewerCanvas component, the coords / number of atoms will not change
const CrystalViewer = ({ frames }: CrystalViewerProps) => {
    const [processedFrames, setProcessedFrames] = useState<ProcessedFrame[]>(
        [],
    );
    const frameNumber = 0; // just have a single frame for now
    useEffect(() => {
        const processedFrames = frames.map((frame): ProcessedFrame => {
            const cornerCoords = getParallelepipedCornerCoords(frame.lattice);
            const middleOfLattice = findMiddleOfCoords(cornerCoords);
            const cartesianCoords = [];
            let furthestDistFromOrigin = 0;
            for (let i = 0; i < frame.fracCoord.length; i++) {
                const fracCoord = frame.fracCoord[i];
                const cartesianCoord = fracCoordToCartesianCoord(
                    fracCoord,
                    frame.lattice,
                );
                const centeredCoords = subtract(
                    cartesianCoord,
                    middleOfLattice,
                );
                furthestDistFromOrigin = Math.max(
                    furthestDistFromOrigin,
                    norm(centeredCoords),
                );
                cartesianCoords.push(centeredCoords);
            }

            // since we shifted all the atoms, we also need to shift the lattice
            const processedLatticePoints = cornerCoords.map((coord) =>
                subtract(coord, middleOfLattice),
            );

            const latticeEdges = getParallelepipedLatticeEdges(
                processedLatticePoints,
            );

            return {
                latticeEdges,
                atomicNumbers: frame.atomicNumbers,
                cartesianCoords: cartesianCoords,
                furthestDistFromOrigin,
            };
        });

        setProcessedFrames(processedFrames);
    }, [frames]);

    if (processedFrames.length === 0) {
        return <div></div>;
    }

    return (
        <CrystalViewerCanvas processedFrame={processedFrames[frameNumber]} />
    );
};

export interface CrystalViewerCanvasProps {
    processedFrame: ProcessedFrame;
}
export interface ProcessedFrame {
    atomicNumbers: number[];
    cartesianCoords: number[][]; // PERF: maybe replace with a fixed size array if there's perf issues
    latticeEdges: Edge[];
    furthestDistFromOrigin: number;
}

const CrystalViewerCanvas = ({ processedFrame }: CrystalViewerCanvasProps) => {
    const mountRef = useRef<HTMLDivElement>(null);

    // I don't think we can cache processed frames and pre-add the atoms them to the scene
    // because if we move the camera and flip between scenes, it'll change?
    // we'll deal with this problem (if it is a caching problem) when we have multiple frames
    useEffect(() => {
        const mount = mountRef.current!;
        const scene = new THREE.Scene();
        // scene.background = new THREE.Color(0xbbbbbb); // Set the background color to white
        scene.background = new THREE.Color(0xffffff); // Set the background color to white
        const camera = new THREE.PerspectiveCamera(
            75,
            mount.clientWidth / mount.clientHeight,
            0.1,
            1000,
        );
        const renderer = new THREE.WebGLRenderer({ antialias: true });
        renderer.setSize(mount.clientWidth, mount.clientHeight);
        // renderer.setSize(window.innerWidth, window.innerHeight);
        renderer.setPixelRatio(window.devicePixelRatio);
        mount.appendChild(renderer.domElement);

        // Add a light to the scene
        const directionalLight = new THREE.DirectionalLight(0xffffff, 1);
        directionalLight.position.set(100, 100, 100).normalize();
        scene.add(directionalLight);

        // Add a second light to the scene so the backside of the atoms is visible
        const directionalLight2 = new THREE.DirectionalLight(0xffffff, 1);
        directionalLight2.position.set(-100, -100, -100).normalize();
        scene.add(directionalLight2);

        // Optionally, add an AmbientLight for softer shadows
        const ambientLight = new THREE.AmbientLight(0xbbbbbb, 0.9); // Soft white light
        scene.add(ambientLight);

        // draw the atoms
        const numAtoms = processedFrame.atomicNumbers.length;
        for (let i = 0; i < numAtoms; i++) {
            const sphere = spheresMap.get(processedFrame.atomicNumbers[i]);
            if (!sphere) {
                console.error(
                    `No sphere found for atomic number ${processedFrame.atomicNumbers[i]}`,
                );
                continue;
            }
            const cartPos = processedFrame.cartesianCoords[i];

            const newSphere = sphere.clone();
            newSphere.position.set(cartPos[0], cartPos[1], cartPos[2]);
            scene.add(newSphere);
        }

        // only draw the edges if the lattice is small enough
        // this is because super large lattices are meant for calculating the energy of molecules, NOT crystal materials
        // the lattice is made large on purpose, so there are no periodic boundary condition interactions
        if (Math.abs(processedFrame.latticeEdges[0][0][0]) < 500) {
            for (const edge of processedFrame.latticeEdges) {
                drawEdge(scene, edge);
            }
        }

        // Camera position
        camera.position.z = 300;

        // OrbitControls for interactivity
        const controls = new OrbitControls(camera, renderer.domElement);
        controls.minDistance = 5; // Set the minimum zoom-in distance
        controls.maxDistance = processedFrame.furthestDistFromOrigin + 8; // Set the maximum zoom-out distance

        // Animation loop
        const animate = () => {
            requestAnimationFrame(animate);
            controls.update();
            renderer.render(scene, camera);
        };

        animate();

        // Clean up on unmount
        return () => {
            mount.removeChild(renderer.domElement);
        };
    }, [processedFrame]);

    return <div ref={mountRef} style={{ width: "100%", height: "100%" }} />;
};

export default CrystalViewer;
