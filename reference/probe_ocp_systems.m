addpath(fullfile(getenv("APPDATA"), "MathWorks", "MATLAB Add-Ons", "Toolboxes", ...
    "ELARA Robot Simulation and Control Toolbox", "examples", "exampleSystems"));

function dumpFrameData(links, outPath)
    MBSim = MBSimulation(links, "displayInfo", false);
    FD = MBSim.MBSys.frameData;
    out = struct();
    out.nDoF = MBSim.MBSys.nDoF;
    out.nFrames = numel(FD);
    out.g0 = MBSim.MBSys.g0;
    out.cSys = MBSim.MBSys.cSys(:)';
    out.dSys = MBSim.MBSys.dSys(:)';
    out.qRef = MBSim.MBSys.qRef(:)';
    out.gConst = MBSim.simPars.g;
    out.nInputs = MBSim.MBSys.nInputs;
    out.indexTCPFrame = MBSim.MBSys.indexTCPFrame;
    out.g_B_TCP = MBSim.MBSys.g_B_TCP;

    fnames = properties(FD);
    fr = struct();
    for k = 1:numel(fnames)
        name = fnames{k};
        try
            val = FD.(name);
            if isnumeric(val) || islogical(val)
                fr.(name) = val;
            end
        catch
        end
    end
    out.frames{1} = fr;

    fid = fopen(outPath, "w");
    fprintf(fid, "%s", jsonencode(out, "PrettyPrint", true));
    fclose(fid);
    fprintf("Wrote %s (nDoF=%d, nFrames=%d, nInputs=%d)\n", outPath, out.nDoF, out.nFrames, out.nInputs);
end

% Planar 2-link pendulum (as used in optimal_control_planar_robot.m)
linksPlanar = systemDef_planarNLinkPendulum("nLinks", 2, "d", 0);
dumpFrameData(linksPlanar, fullfile(pwd, "ocp_planar_frame_data.json"));

% Rigid 3-link lab robot (as used in optimal_control_rigid_robot.m)
linksRigid = systemDef_rigid_robot("d", 0);
dumpFrameData(linksRigid, fullfile(pwd, "ocp_rigid_frame_data.json"));

disp("DONE.");
