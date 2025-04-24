import os
import subprocess

import pymol2

# Dictionary of structure pairs
structure_pairs = {
    "bcl2": (
        "./benchmarking_data/controls/bcl2/viral/EF-CAD53396.1_127_relaxed.pdb",
        "./benchmarking_data/human_structures/AF-P10415-F1-model_v4.pdb",
    ),
    "c1lpt1": (
        "./benchmarking_data/controls/c1lpt1/viral/CF-AAO89306.1_12013_relaxed_pt1.pdb",
        "./benchmarking_data/human_structures/AF-Q8WXC3-F1-model_v4.pdb",
    ),
    "c1lpt2": (
        "./benchmarking_data/controls/c1lpt2/viral/CF-AAO89306.1_12013_relaxed_pt2.pdb",
        "./benchmarking_data/human_structures/AF-P10415-F1-model_v4.pdb",
    ),
    "c4bp": (
        "./benchmarking_data/controls/c4bp/viral/EF-AAO89304.1_12013_relaxed.pdb",
        "./benchmarking_data/human_structures/AF-P04003-F1-model_v4.pdb",
    ),
    "ccr1": (
        "./benchmarking_data/controls/ccr1/viral/EF-AAR31716.1_100_relaxed.pdb",
        "./benchmarking_data/human_structures/AF-P32246-F1-model_v4.pdb",
    ),
    "ccr2": (
        "./benchmarking_data/controls/ccr2/viral/EF-AAD46503.1_158_relaxed.pdb",
        "./benchmarking_data/human_structures/AF-P41597-F1-model_v4.pdb",
    ),
    "cd47": (
        "./benchmarking_data/controls/cd47/viral/CF-AAO89441.1_12013_relaxed.pdb",
        "./benchmarking_data/human_structures/AF-Q08722-F1-model_v4.pdb",
    ),
    "chemokine": (
        "./benchmarking_data/controls/chemokine/viral/EF-AAC55276.1_12002_relaxed.pdb",
        "./benchmarking_data/human_structures/AF-Q9NRJ3-F1-model_v4.pdb",
    ),
    "eif2a": (
        "./benchmarking_data/controls/eif2a/viral/CF-AAO89313.1_12013_relaxed_vaccinia.pdb",
        "./benchmarking_data/human_structures/AF-P05198-F1-model_v4.pdb",
    ),
    "helicase": (
        "./benchmarking_data/controls/helicase/viral/CF-AAO42519.1.1.4_7222_relaxed.pdb",
        "./benchmarking_data/human_structures/AF-Q08211-F1-model_v4.pdb",
    ),
    "ifngr": (
        "./benchmarking_data/controls/ifngr/viral/CF-AAO89469.1_12013_relaxed.pdb",
        "./benchmarking_data/human_structures/AF-P15260-F1-model_v4.pdb",
    ),
    "il10": (
        "./benchmarking_data/controls/il10/viral/CF-CAD53385.1_127_relaxed_ebv.pdb",
        "./benchmarking_data/human_structures/AF-P22301-F1-model_v4.pdb",
    ),
    "il18bp": (
        "./benchmarking_data/controls/il18bp/viral/EF-AAC55182.1_12002_relaxed.pdb",
        "./benchmarking_data/human_structures/AF-O95998-F1-model_v4.pdb",
    ),
    "kinase": (
        "./benchmarking_data/controls/kinase/viral/CF-CAD53438.2_127_relaxed.pdb",
        "./benchmarking_data/human_structures/AF-Q00535-F1-model_v4.pdb",
    ),
    "lfg4": (
        "./benchmarking_data/controls/lfg4/viral/EF-AAL73713.1_12006_relaxed.pdb",
        "./benchmarking_data/human_structures/AF-Q9HC24-F1-model_v4.pdb",
    ),
    "nsp16": (
        "./benchmarking_data/controls/nsp16/viral/CF-QHD43415.1.1.36_10195_relaxed.pdb",
        "./benchmarking_data/human_structures/AF-Q9UI43-F1-model_v4.pdb",
    ),
    "nsp5": (
        "./benchmarking_data/controls/nsp5/viral/CF-QHD43415.1.1.21_10195_relaxed.pdb",
        "./benchmarking_data/human_structures/AF-Q2L4Q9-F1-model_v4.pdb",
    ),
}

# Define custom colors
custom_colors = {"seaweed": [0.314, 0.533, 0.773], "rose": [0.9725, 0.5961, 0.6824]}

# Toggle for high-quality rendering
high_quality = True  # Set to True for ray-traced final renders

# Start PyMOL session
with pymol2.PyMOL() as pymol:
    cmd = pymol.cmd

    for pair_name, (file1, file2) in structure_pairs.items():
        print(f"Processing pair: {pair_name}")

        # Load structures
        mol1 = f"{pair_name}_mol1"
        mol2 = f"{pair_name}_mol2"
        cmd.load(file1, mol1)
        cmd.load(file2, mol2)

        # Align second structure to first
        cmd.cealign(mol1, mol2)

        # Define custom colors
        for name, rgb in custom_colors.items():
            cmd.set_color(name, rgb)

        # Apply visual styles
        cmd.hide("everything")
        cmd.show("cartoon")
        cmd.set("bg_rgb", [1, 1, 1])
        cmd.set("ray_opaque_background", 1)
        cmd.set("cartoon_side_chain_helper", 1)
        cmd.set("ray_shadow", 0)
        cmd.set("stick_radius", 0.15)
        cmd.set("cartoon_loop_cap", 2)
        cmd.set("cartoon_smooth_loops", 1)
        cmd.set("cartoon_loop_quality", 10)
        cmd.set("cartoon_loop_radius", 0.16)
        cmd.set("cartoon_rect_width", 0.3)
        cmd.set("cartoon_fancy_helices", 1)
        cmd.set("cartoon_fancy_sheets", 1)
        cmd.set("cartoon_tube_radius", 0.4)
        cmd.set("cartoon_ring_mode", 3)
        cmd.set("cartoon_ring_width", 0.1)
        cmd.set("cartoon_ring_finder", 1)
        cmd.set("cartoon_ladder_mode", 1)
        cmd.set("cartoon_nucleic_acid_mode", 4)
        cmd.set("cartoon_ring_transparency", 0.5)
        # cmd.set("ray_trace_mode", 1)
        cmd.set("antialias", 1)

        # Apply coloring
        cmd.color("rose", mol1)
        cmd.color("seaweed", mol2)

        # Center and orient view
        cmd.orient()

        # Setup movie frames manually with rotation
        cmd.mset("1 x90")
        for i in range(1, 91):
            cmd.turn("y", 4)
            cmd.mview("store", i)
        cmd.mview("interpolate")
        print("Stored and interpolated 90 frames.")
        print("Frame count reported by PyMOL:", cmd.count_frames())

        # Prepare output folder
        frame_folder = f"{pair_name}_frames_cealign"
        os.makedirs(frame_folder, exist_ok=True)

        # Export PNG frames
        print("Starting PNG frame export...")
        for i in range(1, 91):
            try:
                print(f"Rendering frame {i}")
                cmd.frame(i)
                cmd.refresh()
                out_file = os.path.join(frame_folder, f"{pair_name}_frame{i:04d}.png")
                if high_quality:
                    cmd.png(out_file, width=2100, height=1900, dpi=300, ray=1)
                else:
                    cmd.png(out_file, width=2100, height=1900, dpi=150)
            except Exception as e:
                print(f"Error at frame {i}: {e}")
                break
                
        # Create MP4 using ffmpeg
        mp4_name = f"{pair_name}_cealignment.mp4"
        subprocess.run(
            [
                "ffmpeg",
                "-framerate",
                "15",
                "-i",
                f"{frame_folder}/{pair_name}_frame%04d.png",
                "-vf",
                "crop=1800:800:0:150", #crops width:height:x_offset:y_offset
                "-c:v",
                "libx264",
                "-pix_fmt",
                "yuv420p",
                "-y",
                mp4_name,
            ]
        )

        cmd.reinitialize()
