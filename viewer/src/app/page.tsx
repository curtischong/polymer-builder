"use client";
import ForceLineChart from "@/components/relax/ForcesGraph";
import {
    Frame,
    PolymerVideoViewer,
} from "@/components/relax/PolymerVideoViewer";
import { MantineProvider } from "@mantine/core";
import { Dropzone, FileWithPath } from "@mantine/dropzone";
import React from "react";
import { toast, ToastContainer } from "react-toastify";
import 'react-toastify/dist/ReactToastify.css';
import "@mantine/core/styles.css";

const MAX_FILE_SIZE_BYTES = 1e9;

export default function RelaxPage(){
    const [frames, setFrames] = React.useState<Frame[]>([]);
    // const [forces, setForces] = React.useState<number[][]>([]);

    const onDrop = React.useCallback((files: FileWithPath[]) => {
        if (files.length > 1) {
            toast.error("Please only upload one file at a time");
            return;
        }
        const file = files[0];
        const reader = new FileReader();
        reader.onload = function (e) {
            const contents = e.target?.result;
            if (!contents) {
                toast.error("Failed to read file");
                return;
            }
            const relaxation = JSON.parse(contents as string);
            const fileFrames = [];
            for (const frame of relaxation.frames) {
                fileFrames.push({
                    atomicNumbers: frame.atomic_nums,
                    coords: frame.coords,
                    forces: frame.forces,
                });
            }
            setFrames(fileFrames);
            // setForces(fileForces);
        };
        reader.readAsText(file);
    }, []);

    return (
      <MantineProvider defaultColorScheme="light">
        <div className="mx-4 my-10">
            {frames.length === 0 && (
              <>
                <h1 className="text-2xl font-bold">Drop your relaxation (or click to select it)</h1>
                <Dropzone
                    className="border-dashed border-slate-400 border-[1px] rounded-xl cursor-pointer pb-10 pt-8 h-[300px]"
                    // onReject={handleReject}
                    maxSize={MAX_FILE_SIZE_BYTES}
                    // if we enable this, we'll get logs of warnings since these extensions are not standard mime types
                    accept={{ "application/json": [".json"] }}
                    // {...props}
                    onDrop={onDrop}
                ></Dropzone>
              </>
            )}
            {/* <p>Note: if it gets laggy, you need to refresh the page and try again. Sorry! I an out of time to fix this</p> */}
            <div className="h-[500px] w-[1000px] m-auto">
              {frames.length > 0 && (<p className="mb-4">frames: {frames.length}</p>)}
                <PolymerVideoViewer frames={frames} />
            </div>
            <div className="mt-32">
                {/* {frames.length > 0 && (
                    <div className="flex flex-col my-4">
                        <div>abs force Z direction</div>
                        <div>
                            {frames[frames.length - 1].forces && (
                                <div>
                                    {frames[frames.length - 1].forces
                                        ?.map((f) => Math.abs(f[2]))
                                        .reduce((a, b) => a + b, 0)}
                                </div>
                            )}
                        </div>
                    </div>
                )} */}
                {frames.length > 0 && frames[0].forces && (
                  <ForceLineChart forceSets={frames.map((f) => f.forces || [])} />
                )}
                {/* {frames.length > 0 &&
                    frames[frames.length - 1].forces?.map((f, i) => (
                        <p key={i}>{f.join(", ")}</p>
                    ))} */}
            </div>
            <ToastContainer />
        </div>
        </MantineProvider>
    );
};
