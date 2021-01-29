% This script is used to register TMA cores coordinates

current_path = pwd;



%%

for row = 'A':'L'
    for col = 1:13
        %coreName = 'I-3';
        
        coreName = sprintf('%c-%i', row, round(col));
                % read two reference images
        HE = imread(sprintf('%s/Core_images/%s/H&E.jpg', current_path, coreName));
        CK56 = imread(sprintf('%s/Core_images/%s/CK56.jpg', current_path, coreName));
        Ki67 = imread(sprintf('%s/Core_images/%s/Ki67.jpg', current_path, coreName));
        
        centroids_folder = dir(sprintf('%s/Tissue_Boundary_Rescaled/%s', current_path, coreName));
        filenumber = size(centroids_folder,1);
        
        % create a folder to save all registered files
        
        % comment the two lines if directories are already exist
        reg_dirPath = sprintf('%s/Tissue_Boundary_Registered/%s', pwd, coreName);
        mkdir(reg_dirPath)
        
        %---------------------- Registration begin -----------------------%
        if(filenumber == 12)
            % register images
            
            for flag = 0:9
                switch flag
                    
                    case 0
                        coreType = 'HE';
                        coords_toReg = load(sprintf('%s/Tissue_Boundary_Rescaled/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;
                    
                        % point registration
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'coords_toReg', '-v6');
                        Registered_Pts = [];
                        
                    case 1
                        coreType = 'CK5-6';
                        coords_toReg = load(sprintf('%s/Tissue_Boundary_Rescaled/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;
                    
                        % point registration
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'coords_toReg', '-v6');
                        Registered_Pts = [];
                    %--------------------- H&E section -------------------%
                    case 2
                        coreType = 'Ki67';
                                                
                        coords_toReg = load(sprintf('%s/Tissue_Boundary_Rescaled/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;
                    
                        % point registration
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'coords_toReg', '-v6');
                        Registered_Pts = [];

                
                    case 3
                    
                        coreType = 'P53';
                                            
                        coords_toReg = load(sprintf('%s/Tissue_Boundary_Rescaled/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;
                        % read image and obtain the transformation matrix
                        moveImg = imread(sprintf('%s/Core_images/%s/%s.jpg', current_path, coreName, coreType));
                        movingReg = rigid_transmatrix(moveImg, HE);
                    
                        % point registration
                    
                        [Registered_Pts(:,1),Registered_Pts(:,2)] = transformPointsForward(movingReg.Transformation,coords_toReg(:,1),coords_toReg(:,2));
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'Registered_Pts', '-v6');
                        
                        Registered_Pts = [];

                    case 4
                        coreType = 'Cyclin';
                                        
                        coords_toReg = load(sprintf('%s/Tissue_Boundary_Rescaled/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;
                        
                        % read image and obtain the transformation matrix
                        moveImg = imread(sprintf('%s/Core_images/%s/%s.jpg', current_path, coreName, coreType));
                        movingReg = rigid_transmatrix(moveImg, Ki67);
                        
                        % point registration
                        [Registered_Pts(:,1),Registered_Pts(:,2)] = transformPointsForward(movingReg.Transformation,coords_toReg(:,1),coords_toReg(:,2));
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'Registered_Pts', '-v6');

                        Registered_Pts = [];

                    case 5
                        coreType = 'P16';
                                   
                        coords_toReg = load(sprintf('%s/Tissue_Boundary_Rescaled/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;

                        % read image and obtain the transformation matrix
                        moveImg = imread(sprintf('%s/Core_images/%s/%s.jpg', current_path, coreName, coreType));
                        movingReg = rigid_transmatrix(moveImg, HE);
                        
                        % point registration
                        [Registered_Pts(:,1),Registered_Pts(:,2)] = transformPointsForward(movingReg.Transformation,coords_toReg(:,1),coords_toReg(:,2));
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'Registered_Pts', '-v6');

                        Registered_Pts = [];
                    %--------------------- CK5-6 section -----------------%
                    case 6
                        coreType = 'P63';

                        coords_toReg = load(sprintf('%s/Tissue_Boundary_Rescaled/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;
                        % read image and obtain the transformation matrix
                        moveImg = imread(sprintf('%s/Core_images/%s/%s.jpg', current_path, coreName, coreType));                    
                        movingReg = rigid_transmatrix(moveImg, CK56);
                    
                        % point registration
                        [Registered_Pts(:,1),Registered_Pts(:,2)] = transformPointsForward(movingReg.Transformation,coords_toReg(:,1),coords_toReg(:,2));
                        
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'Registered_Pts', '-v6');

                        Registered_Pts = [];
                        
                    case 7
                    
                        coreType = 'CK20';

                        coords_toReg = load(sprintf('%s/Tissue_Boundary_Rescaled/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;

                        % read image and obtain the transformation matrix
                        moveImg = imread(sprintf('%s/Core_images/%s/%s.jpg', current_path, coreName, coreType));
                        movingReg = rigid_transmatrix(moveImg, CK56);
                    
                        % point registration
                    
                        [Registered_Pts(:,1),Registered_Pts(:,2)] = transformPointsForward(movingReg.Transformation,coords_toReg(:,1),coords_toReg(:,2));
                       
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'Registered_Pts', '-v6');

                        Registered_Pts = [];
                        
                    case 8
                        coreType = 'Her2Neu';

                        coords_toReg = load(sprintf('%s/Tissue_Boundary_Rescaled/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;
   
                        % read image and obtain the transformation matrix
                        moveImg = imread(sprintf('%s/Core_images/%s/%s.jpg', current_path, coreName, coreType));
                        movingReg = rigid_transmatrix(moveImg, CK56);
                        
                        % point registration
                        [Registered_Pts(:,1),Registered_Pts(:,2)] = transformPointsForward(movingReg.Transformation,coords_toReg(:,1),coords_toReg(:,2));

                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'Registered_Pts', '-v6');
                        
                        Registered_Pts = [];
                            
                    case 9
                        coreType = 'GATA3';

                        coords_toReg = load(sprintf('%s/Tissue_Boundary_Rescaled/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;
                  
                        % read image and obtain the transformation matrix
                        moveImg = imread(sprintf('%s/Core_images/%s/%s.jpg', current_path, coreName, coreType));
                        movingReg = rigid_transmatrix(moveImg, CK56);
                        
                        % point registration
                        [Registered_Pts(:,1),Registered_Pts(:,2)] = transformPointsForward(movingReg.Transformation,coords_toReg(:,1),coords_toReg(:,2));
                         
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'Registered_Pts', '-v6');

                        Registered_Pts = [];
                end
                
            end
            
                
        end
            

    end
end

