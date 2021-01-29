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
        
        centroids_folder = dir(sprintf('%s/IHC_Rescaled_Coords/%s', current_path, coreName));
        filenumber = size(centroids_folder,1);
        
        % create a folder to save all registered files
        
        % comment the two lines if directories are already exist
        reg_dirPath = sprintf('%s/IHC_Registered_Coords/%s', pwd, coreName);
        mkdir(reg_dirPath)
        
        %---------------------- Registration begin -----------------------%
        if(filenumber == 11)
            % register images
            
            for flag = 0:8
                switch flag
                    case 0
                        coreType = 'CK56';
                        coords_toReg = load(sprintf('%s/IHC_Rescaled_Coords/%s/CK5-6.mat', current_path, coreName)).coords_toReg;
                    
                        % point registration
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'coords_toReg', '-v6');
                        Registered_Pts = [];
                    %--------------------- H&E section -------------------%
                    case 1
                        coreType = 'Ki67';
                                                
                        coords_toReg = load(sprintf('%s/IHC_Rescaled_Coords/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;
                        % read image and obtain the transformation matrix
                        moveImg = imread(sprintf('%s/Core_images/%s/%s.jpg', current_path, coreName, coreType));                    
                        movingReg = rigid_transmatrix(moveImg, HE);
                    
                        % point registration
                        [Registered_Pts(:,1),Registered_Pts(:,2)] = transformPointsForward(movingReg.Transformation,coords_toReg(:,1),coords_toReg(:,2));
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'Registered_Pts', '-v6');

                        Registered_Pts = [];

                
                    case 2
                    
                        coreType = 'P53';
                                            
                        coords_toReg = load(sprintf('%s/IHC_Rescaled_Coords/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;
                        % read image and obtain the transformation matrix
                        moveImg = imread(sprintf('%s/Core_images/%s/%s.jpg', current_path, coreName, coreType));
                        movingReg = rigid_transmatrix(moveImg, HE);
                    
                        % point registration
                    
                        [Registered_Pts(:,1),Registered_Pts(:,2)] = transformPointsForward(movingReg.Transformation,coords_toReg(:,1),coords_toReg(:,2));
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'Registered_Pts', '-v6');
                        
                        Registered_Pts = [];

                    case 3
                        coreType = 'Cyclin';
                                        
                        coords_toReg = load(sprintf('%s/IHC_Rescaled_Coords/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;
                        
                        % read image and obtain the transformation matrix
                        moveImg = imread(sprintf('%s/Core_images/%s/%s.jpg', current_path, coreName, coreType));
                        movingReg = rigid_transmatrix(moveImg, HE);
                        
                        % point registration
                        [Registered_Pts(:,1),Registered_Pts(:,2)] = transformPointsForward(movingReg.Transformation,coords_toReg(:,1),coords_toReg(:,2));
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'Registered_Pts', '-v6');

                        Registered_Pts = [];

                    case 4
                        coreType = 'P16';
                                   
                        coords_toReg = load(sprintf('%s/IHC_Rescaled_Coords/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;

                        % read image and obtain the transformation matrix
                        moveImg = imread(sprintf('%s/Core_images/%s/%s.jpg', current_path, coreName, coreType));
                        movingReg = rigid_transmatrix(moveImg, HE);
                        
                        % point registration
                        [Registered_Pts(:,1),Registered_Pts(:,2)] = transformPointsForward(movingReg.Transformation,coords_toReg(:,1),coords_toReg(:,2));
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'Registered_Pts', '-v6');

                        Registered_Pts = [];
                    %--------------------- CK5-6 section -----------------%
                    case 5
                        coreType = 'P63';

                        coords_toReg = load(sprintf('%s/IHC_Rescaled_Coords/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;
                        % read image and obtain the transformation matrix
                        moveImg = imread(sprintf('%s/Core_images/%s/%s.jpg', current_path, coreName, coreType));                    
                        movingReg = rigid_transmatrix(moveImg, CK56);
                    
                        % point registration
                        [Registered_Pts(:,1),Registered_Pts(:,2)] = transformPointsForward(movingReg.Transformation,coords_toReg(:,1),coords_toReg(:,2));
                        
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'Registered_Pts', '-v6');

                        Registered_Pts = [];
                        
                    case 6
                    
                        coreType = 'CK20';

                        coords_toReg = load(sprintf('%s/IHC_Rescaled_Coords/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;

                        % read image and obtain the transformation matrix
                        moveImg = imread(sprintf('%s/Core_images/%s/%s.jpg', current_path, coreName, coreType));
                        movingReg = rigid_transmatrix(moveImg, CK56);
                    
                        % point registration
                    
                        [Registered_Pts(:,1),Registered_Pts(:,2)] = transformPointsForward(movingReg.Transformation,coords_toReg(:,1),coords_toReg(:,2));
                       
                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'Registered_Pts', '-v6');

                        Registered_Pts = [];
                        
                    case 7
                        coreType = 'Her2Neu';

                        coords_toReg = load(sprintf('%s/IHC_Rescaled_Coords/%s/%s.mat', current_path, coreName, coreType)).coords_toReg;
   
                        % read image and obtain the transformation matrix
                        moveImg = imread(sprintf('%s/Core_images/%s/%s.jpg', current_path, coreName, coreType));
                        movingReg = rigid_transmatrix(moveImg, CK56);
                        
                        % point registration
                        [Registered_Pts(:,1),Registered_Pts(:,2)] = transformPointsForward(movingReg.Transformation,coords_toReg(:,1),coords_toReg(:,2));

                        save(sprintf('%s/%s.mat', reg_dirPath, coreType),'Registered_Pts', '-v6');
                        
                        Registered_Pts = [];
                            
                    case 8
                        coreType = 'GATA3';

                        coords_toReg = load(sprintf('%s/IHC_Rescaled_Coords/%s/uv-GATA3.mat', current_path, coreName)).coords_toReg;
                  
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

