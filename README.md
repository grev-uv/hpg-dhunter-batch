# HPG-Dhunter batch process

If you want to use this tool just now, there is an executable file for Linux x86_64 systems. This compressed file is available at [releases page](../../releases). But, previously, CUDA and Nvidia drivers must be installed in your system. For that, you can go to [System requirements section](README.MD#L81)

**HPG-Dhunter batch**  is a complementary tool of [**HPG-Dhunter**](https://github.com/grev-uv/hpg-dhunter) tool that automatically detects all the Differentially Methylated Regions (DMRs) among the considered samples for all the chromosomes in the genome. It is also based on the Discrete Wavelet Transform (DWT), and it provides a list of all the DMRs found.

**HPG-Dhunter batch** identifier is a powerful tool that uses the high performance parallel computing capabilites of GPUs and the CUDA programming interface model to detect DMRs and save the results into a .gff and .csv file, minimizing the CPU-GPU communication. The batch process allows the identification of DMRs analyzing all the samples together. 

## Handling
**HPG-Dhunter batch** shows a user interface (UI) whose design has been developed according to the usability principles.The DMR detection process follows a pipeline that begins selecting the cases and control files. After that, the ratio between the methylated coverage and the total coverage over each chromosome position is calculated and upload to the global memory of the GPU device in batches. With all the the computed results, it is possible to identify the DMRs among all the selected samples.

This is the UI:

![](images/interface_batch.png)

where:
1. Select the chromosome of reference.
2. Open directory browser to select the directory of each case sample.
3. Delete the selected file.
4. Climb one position the selected file.
5. Low one position the selected file
6. List of selected files to analyze
7. Open directory browser to select the directory of each control sample.
8. Open directory browser to select the directory where save the file with DMRs
9. Select the kind of signal to analyze (mC and/or hmC)
10. Select the number of chromosomes to analyze (all / list).
11. Progress bar.
12. Some system information
13. Button to start the process.
14. Button to stop the process.
15. Slider to select the wavelet transform level
16. Slider to select the threshold for DMR identification
17. Slider to select the minimum ratio of methylated positions by region
18. Slider and text window to select the minimum coverage of hmC signal
19. Slider and text window to select the minimum coverage of mC signal
20. Text window to write de list of chromosomes to analyze.
21. Select the direction to analyze (forward and/or reverse).
22. Select the kind of analyze (grouped or single samples)
23. Text window showing the path to save the results.

There is an important change to do before launching the compilation. The file [hpg_dhunter.pro](src/hpg_dhunter.pro#L58) needs the path to cuda sdk installation at line 58:
```
CUDA_DIR = /path/to/cuda/sdk/cuda
```

## System requirements
The HPG-Dhunter-batch tool, as a complementary tool of HPG-Dhunter visualizer, is the next step after HPG-HMapper tool. Therefore, the system requirements are the same as the ones for HPG-HMapper, plus a GPU device.
HPG-Dhunter should work properly in a station with the following set-up:
- A 64 bit Intel CPU compatible with SSE4.2.
- The DNA data for DMR tasks needs as much adjacent memory as the number of samples by the length of the largest chromosome to be analized. This parameter has a direct relationship with the global memory available in the GPU device. The test was done with 32 MB of RAM.
- The amount of samples that HPG-Dhunter can analize at the same time directly depends on the amount of the device memory. Working with a Nvidia GeForce GTX 1080 with 8 GB of GRAM, it is possible to analyze and visualize up to six samples of human chromosome-21 or up to four human chromosome-10, or up to two human chromosome-1 at the same time.
- The CUDA compilation is configured to a single device with Nvidia Pascal GPU architecture. So, the devices that will work properly are Titan XP and X models, Tesla P40, P6 and P4 models, Quadro P6000, P5000, P4000 models, GeForce GTX 1080Ti, 1080, 1070Ti, 1070 models, and others easy to find here.
- The Nvidia driver is needed (v384 or higher).
- The CUDA API is needed(v9 or higher).
- For a complete installation of Nvidia drivers and CUDA, it can download just one file from [here](https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64) and follow the instructions, that appear in the same page after select the distribution and version, for a pretty installation.

## Build
The way to build HPG-Dhunter batch identifier in your system is opening the software as a project inside an installed QtCreator (> v4.5, Qt > v5.8, GCC 5) IDE and build it from there.

In the next future, another available way will be to handling this software as a cloud service.

## Issues
HPG-Dhunter - Copyright (C) 2018 - grev-uv
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it under certain conditions; visit https://www.gnu.org/copyleft/gpl.html for details.
However, if you find any bugs, issues, want a specific feature added or need help, feel free to add an issue or extend an existing one. Pull requests are welcome.


## License
HPG-Dhunter batch identifier is free software and licensed under the GNU General Public License version 3.
HPG-Dhunter batch identifier was developed under Qt as a platform application development framework for linux/ubuntu desktop, using a free software LGPL v3 license.

## Contact
Contact any of the following developers for any enquiry:
- Juanma Orduña (juan.orduna@uv.es). 
- Mariano Pérez (mariano.perez@uv.es). 
- Lisardo Fernández (lisardo.fernandez@uv.es). 
