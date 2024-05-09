
%-------------------------------------------------------------------------
                             % Main Function %
%-------------------------------------------------------------------------

function student_number = tp1_122422_rumore()
    
    % In this function is contained the code which was used to 
    % developed the TP1 project. All the code is comment to have
    % a better understanding of all the development process. 

    %First of all we add the path of the previous folder

    addpath('../');

    %To read all the images i call the dir function
    
    images_in_previous_directory = dir('../svpi2024_TP1_img_*.png');

    %each element of the list that we creates is now a struct. We want a
    %list of images to work with. We are gonna create the list of images
    %thanks to the method create_path

    list_of_all_images = {};

    for i = 1:size(images_in_previous_directory)

        %i take tha path of the current image
        current_path = create_path(images_in_previous_directory(i));
        
        %i read the current image
        current_image = imread(current_path);

        %i add the image to the list
        list_of_all_images{i} = current_image;

    end;

    %Now i'm gonna create the matrix that will contain all of my results
    

    rows = size(list_of_all_images, 1); 
    columns = 17; %number or statistics that are requested

    %I will now associate an index number to each statistic to be able to
    %remeber easily at the columns of the matrix. 

    % NumMec = 1;
    % NumSeq = 2;
    % NumImg = 3;
    % TotNm = 4;
    % TotCB = 5;
    % TotQR = 6;
    % R0 = 7;
    % R90 = 8;
    % R180 = 9;
    % R270 = 10;
    % ReflCB = 11;
    % BadCB = 12;
    % TotDigCB = 13;
    % CBL = 14;
    % CBR = 15;
    % CBG = 16;
    % StringCB = 17;

    statistics_matrix = zeros(rows, columns - 1);

    %for the stirng i need another structure

    central_digits(60, 1) = struct();

    for i = 1:60
        central_digits(i).string_number = "";
    end

    %I already know that in the first column i have to put my student
    %number so i'm gonna do it right now
    %same thing for NumSeq and NumImg
    
    for i = 1:size(images_in_previous_directory)

        %adding the NumMec
        statistics_matrix(i, 1) = 122422;

        %adding NumSeq
        if(i <= 20)
            statistics_matrix(i, 2) = 330;
        elseif (i > 20 & i <= 40)
            statistics_matrix(i, 2) = 420;
        else
            statistics_matrix(i, 2) = 530;
        end
            
        %adding NumImg
        statistics_matrix(i, 3) = mod(i-1, 20) + 1;
    end

    %Now for all the images we call a function which will handle all the
    %images

    for i = 1:size(list_of_all_images, 2)
       
       [statistics_matrix, central_digits] = start_processing(list_of_all_images{i}, statistics_matrix, central_digits, i);

    end

    %I'm gonna sort all the strings that i obtained from de processing

    for i = 1:size(list_of_all_images, 2)
        
        array_of_char = char(central_digits(i).string_number);
        sorted_array = sort(array_of_char);
        central_digits(i).string_number = string(sorted_array);
        central_digits(i).string_number

    end

    %creating the txt file

    file_of_results = fopen('tp1_122422.txt', 'w');

    for i = 1:size(statistics_matrix, 1)

        result_row = ' ';

        for j = 1:size(statistics_matrix, 2)

            result_row = [result_row, num2str(statistics_matrix(i, j))];

            if j < size(statistics_matrix, 2)

                result_row = [result_row, ',']; 

            end

        end

        struct_element = central_digits(i).string_number;

        results_row = strcat(result_row, ',', struct_element);

        fprintf(file_of_results, '%s\n', results_row);
    end

    fclose(file_of_results)
    
    statistics_matrix

    % Only the student number is returned as requested.
    student_number = 122422;

end



%-------------------------------------------------------------------------
                             % Auxiliar Functions %
%-------------------------------------------------------------------------


function path_of_the_image = create_path(info_image)

    %This method takes in input the image infos to create the path
    %in order to save the image
    
    image_name = info_image.name;
    image_folder = info_image.folder;
    path_of_the_image = fullfile(image_folder,'/', image_name);

end



function [new_statistics_matrix, new_central_digits] = start_processing(image, statistics_matrix, central_digits, index)

    %In this function we will start the process of the images. This will be
    %the scheme:

    % 1. Finding the sub images
    % 2. Filter all the sub images
    % 3. Check if a sub images is a bar code
    %   3.1 If is a bar code decode the numbers
    %   3.2 If is not  bar code we check if it is a QR code
    %       3.2.1 If is not a bar or a qr code is a object with no meaning
    
    
    %1. Finding the sub images
    sub_images2 = vs_getsubimages(image);
    sub_images = segmentation(image);

    %2. Filter the sub images
    sub_images_filtered = {};

   for i = 1:size(sub_images, 2)
        
        new_image = filtered_image(sub_images{i});
     
        sub_images_filtered{i} = new_image;
  
    end

    for j = 1:size(sub_images_filtered, 2)
        
        %3.1 Veirify if is a bar code
        [valid_object, statistics_matrix, central_digits] = bar_code_verify(sub_images_filtered{j}, statistics_matrix, central_digits, index);

        %In the previous method we call also the method to verify if it is
        %a qr code. 

        %If we obtained a non valid object we update the number of TotNm = 4;
        if valid_object == false

            statistics_matrix(index, 4) = statistics_matrix(index, 4) + 1;

        end
    end 

    %we update the matrix
    new_central_digits = central_digits;
    new_statistics_matrix = statistics_matrix;

end

function array_of_sub_images = segmentation(image)
    % This function segments the input image into its large components and returns each as a sub-image, respecting a minimum area of 20,000 pixels for each component

    
    % Create an image where non-border pixels (not 0 or 255) are set to maximum intensity (white)
    no_gray_image = 255 * ones(size(image), 'like', image);
    no_gray_image(image == 0 | image == 255) = image(image == 0 | image == 255);

    % Step 2: Edge detection
    % Using the Canny edge detector
    edges = edge(no_gray_image, 'Sobel');

    % Step 3: Find connected components from edge image
    % Inverse the edge image to get the connected regions
    filled = imfill(edges, 'holes');
    imshow(filled);

    % Label connected components
    [labeledImage, numberOfBlobs] = bwlabel(filled, 8);

    

    % Step 4: Extract and store each connected component as a sub-image if it meets the area requirement
    % Preallocate a cell array to hold the sub-images
    array_of_sub_images = {};
    component_count = 0;



% Passo aggiuntivo: Calcolo della lunghezza del perimetro
% e filtraggio delle componenti basato sulla lunghezza
for k = 1:numberOfBlobs
    % Crea un'immagine binaria per ciascuna componente
    componentImage = labeledImage == k;

    % Calcola il perimetro della componente
    perimeter = bwperim(componentImage);
    perimeterLength = sum(perimeter(:));

    % Controlla l'area e la lunghezza del perimetro della componente
    componentArea = sum(componentImage(:));
    if componentArea >= 2000 && componentArea <= 250000 && ...
       perimeterLength >= 40 && perimeterLength <= 20000
        % Se la componente soddisfa i requisiti di area e lunghezza del perimetro
        component_count = component_count + 1;
        % Usa regionprops per ottenere il bounding box di ciascuna componente
        stats = regionprops(componentImage, 'BoundingBox');

        % Ritaglia la componente dall'immagine originale
        boundingBox = ceil(stats.BoundingBox);
        x1 = boundingBox(1);
        y1 = boundingBox(2);
        x2 = x1 + boundingBox(3) - 1;
        y2 = y1 + boundingBox(4) - 1;
        subImage = image(y1:y2, x1:x2);

        % Salva la sub-immagine nell'array di celle
        array_of_sub_images{component_count} = subImage;
    end
end
    % Optionally, display each sub-image
    % for idx = 1:component_count
    %     imshow(array_of_sub_images{idx}), title(sprintf('Component %d', idx));
    % end

end


function output_image = filtered_image(input_image)

    %I will not consider the border of the image
    cropped_image = input_image(3:end-2, 3:end-2);

    % this operation is useful to understand if we need to do a more
    % complex operation on the image while filtering or smìimpky binarize
    % it
    intensity_levels = unique(cropped_image);
    count_levels = numel(intensity_levels);

    if count_levels > 3

        thresholds = multithresh(cropped_image, 2);
        quantized_image = imquantize(cropped_image, thresholds);
        dark_region_mask = quantized_image == 1;
       
        if ndims(cropped_image) == 2
            dark_regions_colored = label2rgb(dark_region_mask);
        else
            dark_regions_colored = cropped_image .* uint8(dark_region_mask);
        end

        gray_image = rgb2gray(dark_regions_colored);
        output_image = imbinarize(gray_image);

    else
        
        contrast_enhanced_image = imadjust(cropped_image);
        output_image = imbinarize(contrast_enhanced_image, 'adaptive', 'ForegroundPolarity', 'dark');

    end

end





function [valid_object, new_statistic_matrix] = qr_code_verify(image, statistic_matrix, index)

    %i initialize the results variables
    valid_object = false;
    new_statistic_matrix = statistic_matrix;

    %first of all we check the dimension, the qr code is a square
    [image_rows, image_columns] = size(image);

    %if it is not a square and it is not a bar code it is an invalid object
    if image_rows ~= image_columns
        
        return;

    end

    %i try to find the quiet zone    
    line_on_top_to_remove = 0;
    for i = 1:image_rows
        if all(image(i, :) == 1)
            line_on_top_to_remove = line_on_top_to_remove + 1;
        else
            break;
        end
    end
    

    column_on_left_to_remove = 0;
    for j = 1:image_columns
        if all(image(:, j) == 1)
            column_on_left_to_remove = column_on_left_to_remove + 1;
        else
            break;
        end
    end

    line_on_bottom_to_remove = 0;
    for i = image_rows:-1:1
        if all(image(i, :) == 1)
            line_on_bottom_to_remove = line_on_bottom_to_remove + 1;
        else
            break;
        end
    end
    
    column_on_right_to_remove = 0;
    for j = image_columns:-1:1
        if all(image(:, j) == 1)
            column_on_right_to_remove = column_on_right_to_remove + 1;
        else
            break;
        end
    end
    
    % we can now remove the quiet zone
    qr_code = image((line_on_top_to_remove + 1):(image_rows - line_on_bottom_to_remove), ...
                    (column_on_left_to_remove+1):(image_columns - column_on_right_to_remove));

    %We can now evaluate the 17+4n rule

    [qr_rows, qr_columns] = size(qr_code);

    for version = 1:40  % controlling the different versions
        expected_modules = 17 + 4 * version;
        if qr_rows == qr_columns && qr_rows == expected_modules
            break;
        end
    end
    
    %if it is a square and it does not respect the 17+4n rule it is an
    %invalid object
    if version > 40

        return;

    end

    %we can now try to find the sub pattern

    qr_code_fill = imfill(qr_code,'holes');
    qr_code_area = bwareaopen(qr_code_fill, 50);
    qr_code_final = bwlabel(qr_code_area);
    statistics = regionprops(qr_code_final,'boundingbox', 'Area');

    %is a qr code if we have at least 4 reagions (including the border)
    %otherwise it is an invalid object
    if numel(statistics) < 4

        return

    end

    %i update the value of qr code
    valid_object = true;
    statistic_matrix(index, 6) = statistic_matrix(index, 6) + 1;
    new_statistic_matrix = statistic_matrix;
    
end


function [valid_object, new_statistics_matrix, new_central_digits] = bar_code_verify(image, statistics_matrix, central_digits, index)

    %i initialize the results variables
    valid_object = false;
    new_central_digits = central_digits;
    new_statistics_matrix = statistics_matrix;

    %i define the start and end sequence to verify if i'm working with a
    %bar code
    start_sequence = [0 0 1 0 1 1 0 1 1 1 0];
    end_sequence = [0 1 1 1 0 0 0 1 0 1 0 0];
    

    [rows, columns] = size(image);
    
    
    %i have to add this control because some of the images without meaning
    %could be total white and this could lead to an error
    mean_white = mean(image(:));

    if mean_white >= 0.98

        new_statistics_matrix = statistics_matrix;
        return;

    end
    
    %i will remove also if the "quiet zone" to work directly with the
    %information zone of the  bar code

    line_on_top_to_remove = 0;
    for i = 1:rows
        if all(image(i, :) == 1)
            line_on_top_to_remove = line_on_top_to_remove + 1;
        else
            break;
        end
    end
    
    column_on_left_to_remove = 0;
    for j = 1:columns
        if all(image(:, j) == 1)
            column_on_left_to_remove = column_on_left_to_remove + 1;
        else
            break;
        end
    end
    

    line_on_bottom_to_remove = 0;
    for i = rows:-1:1
        if all(image(i, :) == 1)
            line_on_bottom_to_remove = line_on_bottom_to_remove + 1;
        else
            break;
        end
    end
    

    column_on_right_to_remove = 0;
    for j = columns:-1:1
        if all(image(:, j) == 1)
            column_on_right_to_remove = column_on_right_to_remove + 1;
        else
            break;
        end
    end
    
    bar_code = image((line_on_top_to_remove + 1):(rows - line_on_bottom_to_remove), ...
                    (column_on_left_to_remove+1):(columns - column_on_right_to_remove));   
    

    % we know that the bat code can have different angles and coul be
    % reflected so we will proceed as follows:

    %1. Try if we can recognize the bar code for one of the angles
    %2. If we did not recognize the bar code we have to make sure it is not
    %reflected. We will try to reflect it and try again for every angle if
    %we can recognize it

    orientations = {'R0', 'R90', 'R180', 'R270', 'R0F', 'R90F', 'R180F', 'R270F'};
    rotations = [0, -90, -180, -270];
    valid_orientation = '';
    found = false;
    pixel_per_bar = 1;

    while ~found && pixel_per_bar <= 3

        %i  update start and end sequence because some bar code could have
        %the bars wider than 1 single pixel
        updated_start_sequence = repelem(start_sequence, pixel_per_bar);
        updated_end_sequence = repelem(end_sequence, pixel_per_bar);

        for i = 1:length(orientations)

            if found

                break;

            end

            for j = 1:length(rotations)

   
                rotated_bar_code = imrotate(bar_code, rotations(j));
                [rows_bar_code_rotated, columns_bar_code_rotated] = size(rotated_bar_code);

                if columns_bar_code_rotated >= (11 * pixel_per_bar) 
                    
                    %i select the first row of the bar code where the
                    %infromation should be. 
                    test_row = rotated_bar_code(1, :);
                    
                    if length(test_row) - (11*pixel_per_bar+1) > 0

                        %i have to change the way i take the start and end
                        %sequence from the test row with respect to the
                        %number of pixel
                        if pixel_per_bar == 1

                            foundStart = isequal(updated_start_sequence, test_row(:, 1:11));
                            foundEnd = isequal(updated_end_sequence, test_row(:, end-11:end));

                        else

                            foundStart = isequal(updated_start_sequence, test_row(:, 1:(11*pixel_per_bar)));
                            foundEnd = isequal(updated_end_sequence, test_row(:, end-(11*pixel_per_bar+1):end));

                        end
            
                        if foundStart && foundEnd

                            valid_orientation = orientations{j};
                      
                            %we have some cases in which the orientation is
                            %found but is not the right one due to
                            %similarity to start and end sequence. For
                            %istance, a R0 bar code have the information
                            %similar to an R180F
                            if (rotations(j) == 0 ) && all(all(image(1:15, 1:15) == 1))

                                break;

                            end

                            %once we found a bar code we have to verify if
                            %it is valid with respect to the R-L-G
                            %codification
                            [invalid, statistics_matrix, central_digits] = decode_number_in_bar_code(rotated_bar_code, pixel_per_bar, statistics_matrix, central_digits, index);
                            found = true;

                            break;
                        end
                    end
                end
            end
        end
    
        %now we try with the bar code flipped
        if ~found
            
            flipped_bar_code = flip(bar_code, 2);


            %we do the same operations as before but with the code flipper
            %on the y axis
            for i = 1:length(orientations)

                if found

                    break;

                end

                for j = 1:length(rotations)
                    
                    rotated_bar_code = imrotate(flipped_bar_code, rotations(j));
                    [rows_bar_code_rotated, columns_bar_code_rotated] = size(rotated_bar_code);
    
                    if columns_bar_code_rotated >= (11 * pixel_per_bar) 

                        test_row = rotated_bar_code(1, :);
                        
                        if length(test_row) - (11*pixel_per_bar+1) > 0

                             if pixel_per_bar == 1

                                foundStart = isequal(updated_start_sequence, test_row(:, 1:11));
                                foundEnd = isequal(updated_end_sequence, test_row(:, end-11:end));

                             else

                                foundStart = isequal(updated_start_sequence, test_row(:, 1:(11*pixel_per_bar)));
                                foundEnd = isequal(updated_end_sequence, test_row(:, end-(11*pixel_per_bar+1):end));

                             end
                
                            if foundStart && foundEnd

                                valid_orientation = orientations{j+4};

                                if (rotations(j) == 0 ) && all(all(image(1:15, 1:15) == 1))

                                    break;

                                end

                                %also in this case we have to be careful
                                %with the orientations
                                if (rotations(j) == -270) && all(all(image(1:5, end-5:end)==1))
                                   
                                    valid_orientation = 'R90F';
                                    found = true;
                                    [invalid, statistics_matrix, central_digits] = decode_number_in_bar_code(rotated_bar_code, pixel_per_bar, statistics_matrix, central_digits, index); 
                                    break;

                                end

                                if (rotations(j) == -90) && all(all(image(1:5, 1:5)==1))

                                    valid_orientation = 'R270F';
                                    found = true;
                                    [invalid, statistics_matrix, central_digits] = decode_number_in_bar_code(rotated_bar_code, pixel_per_bar, statistics_matrix, central_digits, index);
                                    break;

                                end

                                [invalid, statistics_matrix, central_digits] = decode_number_in_bar_code(rotated_bar_code, pixel_per_bar, statistics_matrix, central_digits, index);
                                found = true;
                                break;

                            end
                        end
                    end
                end
            end
        end

        %if even without the flipped image we found a correspondance we
        %have to try to consider more pixel per information bar
        pixel_per_bar = pixel_per_bar +1;

    end

    

    %if we have found some valid orientation we update the matrix
    if ~isempty(valid_orientation)
       
        valid_object = true; 

        if strcmp(valid_orientation, 'R0')

            statistics_matrix(index, 7) = statistics_matrix(index, 7) + 1;
            statistics_matrix(index, 5) = statistics_matrix(index, 5) + 1; 

        elseif strcmp(valid_orientation, 'R90')

            statistics_matrix(index, 8) = statistics_matrix(index, 8) + 1;
            statistics_matrix(index, 5) = statistics_matrix(index, 5) + 1; 

        elseif strcmp(valid_orientation, 'R180')

            statistics_matrix(index, 9) = statistics_matrix(index, 9) + 1;
            statistics_matrix(index, 5) = statistics_matrix(index, 5) + 1;
            
        elseif strcmp(valid_orientation, 'R270')

            statistics_matrix(index, 10) = statistics_matrix(index, 10) + 1;
            statistics_matrix(index, 5) = statistics_matrix(index, 5) + 1;
            
        elseif strcmp(valid_orientation, 'R0F')

            statistics_matrix(index, 7) = statistics_matrix(index, 7) + 1;
            statistics_matrix(index, 5) = statistics_matrix(index, 5) + 1;
            
        elseif strcmp(valid_orientation, 'R90F')

            statistics_matrix(index, 8) = statistics_matrix(index, 8) + 1;
            statistics_matrix(index, 5) = statistics_matrix(index, 5) + 1;
            
        elseif strcmp(valid_orientation, 'R180F')

            statistics_matrix(index, 9) = statistics_matrix(index, 9) + 1;
            statistics_matrix(index, 5) = statistics_matrix(index, 5) + 1;
           
        elseif strcmp(valid_orientation, 'R270F')

            statistics_matrix(index, 10) = statistics_matrix(index, 10) + 1;
            statistics_matrix(index, 5) = statistics_matrix(index, 5) + 1;

        end

        %if the bar code is valid and is reflected we update the 11th
        %position in the matrix
        if ~invalid && (strcmp(valid_orientation, 'R0F') | strcmp(valid_orientation, 'R90F') | strcmp(valid_orientation, 'R180F') | strcmp(valid_orientation, 'R270F'))

            statistics_matrix(index, 11) = statistics_matrix(index, 11) + 1;
       
        end
    
    else

        %if it is not a valid bar_code i chack if it is a valid qr_code
        [valid_object, new_statistics_matrix] = qr_code_verify(image, statistics_matrix, index);
        return

    end
    
    new_statistics_matrix = statistics_matrix;
    new_central_digits = central_digits;
    
end


function [invalid, new_statistic_matrix, new_central_digits] = decode_number_in_bar_code(image, pixel_per_bar, statistic_matrix, central_digits, index)

    %i initialize the results variables
    new_statistic_matrix = statistic_matrix;
    new_central_digits = central_digits;
    invalid = 0;

            L = [
                1 1 1 0 0 1 0; 
                1 1 0 0 1 1 0; 
                1 1 0 1 1 0 0; 
                1 0 0 0 0 1 0;
                1 0 1 1 1 0 0;
                1 0 0 1 1 1 0; 
                1 0 1 0 0 0 0;
                1 0 0 0 1 0 0; 
                1 0 0 1 0 0 0; 
                1 1 1 0 1 0 0;
            ];

            R = [
                0 0 0 1 1 0 1; 
                0 0 1 1 0 0 1; 
                0 0 1 0 0 1 1; 
                0 1 1 1 1 0 1; 
                0 1 0 0 0 1 1; 
                0 1 1 0 0 0 1; 
                0 1 0 1 1 1 1; 
                0 1 1 1 0 1 1; 
                0 1 1 0 1 1 1; 
                0 0 0 1 0 1 1; 
            ];

            G = [
                1 0 1 1 0 0 0; 
                1 0 0 1 1 0 0; 
                1 1 0 0 1 0 0; 
                1 0 1 1 1 1 0; 
                1 1 0 0 0 1 0; 
                1 0 0 0 1 1 0; 
                1 1 1 1 0 1 0; 
                1 1 0 1 1 1 0; 
                1 1 1 0 1 1 0; 
                1 1 0 1 0 0 0; 
            ];

    %we are going to consider only the part in the middle of the bar code
    %because we already verified the initial and final sequence
    
    internal_bar_code = image(:, (11*pixel_per_bar)+1:end-(11*pixel_per_bar) - pixel_per_bar);
    
    %i will compute the number incoded
    numbers = size(internal_bar_code,2) / (7*pixel_per_bar);

    %this string will be used to store the number of the bar code
    string_of_number = "";


    for i = 1:numbers

        %i extract the segment that represents the current digit
        digit_start = (i-1) * 7 *pixel_per_bar + 1;
        digit_end = i * 7 *pixel_per_bar;
        digit_segment = internal_bar_code(:, digit_start:digit_end);
        digit_segment = digit_segment(1,:);
        
        %i will use boolean values to store the current codification
        is_L = false;
        is_R = false;
        is_G = false;

        
        %now i simply loop on every codification and see if a match exists
        for j = 1:10 

            l_code = L(j, :);
            r_code = R(j, :);
            g_code = G(j, :);

            l_repeated = repelem(l_code, pixel_per_bar);
            r_repeated = repelem(r_code, pixel_per_bar);
            g_repeated = repelem(g_code, pixel_per_bar);

            if isequal(digit_segment, l_repeated )

                is_L = true;

                if i == 1
                    
                    %the first time i encounter a L codification i update
                    %the matrix
                    statistic_matrix(index, 14) = statistic_matrix(index, 14) +1;

                end

                %storing relevant information
                statistic_matrix(index, 13) = statistic_matrix(index, 13) + 1;
                string_of_number = strcat(string_of_number, num2str(j-1));
                break; 

            end
            

            if isequal(digit_segment, r_repeated)

                is_R = true;

                if i == 1

                    statistic_matrix(index, 15) = statistic_matrix(index, 15) +1;

                end

                statistic_matrix(index, 13) = statistic_matrix(index, 13) + 1;
                string_of_number = strcat(string_of_number, num2str(j-1));
                break; 

            end
            

            if isequal(digit_segment, g_repeated)

                is_G = true;

                if i == 1

                    statistic_matrix(index, 16) = statistic_matrix(index, 16) +1;

                end

                statistic_matrix(index, 13) = statistic_matrix(index, 13) + 1;
                string_of_number = strcat(string_of_number, num2str(j-1));
                break; 

            end
        end
                
        % If the segment doesn't match any codification, the barcode is invalid
        if ~(is_L || is_R || is_G) 

            invalid = 1;
            statistic_matrix(index, 12) = statistic_matrix(index, 12) +1;
            new_statistic_matrix = statistic_matrix;
            return; 

        end

    end

    %i extract the central digit of the sequence
    position_central_digit = floor(strlength(string_of_number)/2)+1;
    central_number = extractBetween(string_of_number, position_central_digit, position_central_digit);
    central_digits(index).string_number = strcat(central_digits(index).string_number, central_number);
    new_central_digits = central_digits;
    new_statistic_matrix = statistic_matrix;

end
