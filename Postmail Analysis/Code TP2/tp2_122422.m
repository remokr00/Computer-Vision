
%-------------------------------------------------------------------------
                             % Main Function %
%-------------------------------------------------------------------------

function student_number = tp2_122422()
 
    % In this function is contained the code which was used to 
    % developed the TP1 project. All the code is comment to have
    % a better understanding of all the development process. 

    
    %First of all we add the path of the previous folder

    addpath('../');

    %loading the neural network
   


    %To read all the images i call the dir function
    
    images_in_previous_directory = dir('../svpi2024_TP2_img_*.png');

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

    end

    %Now i'm gonna create the matrix that will contain all of my results
    

    rows = size(list_of_all_images, 1); 
    columns = 15; %number or statistics that are requested

    %I will now associate an index number to each statistic to be able to
    %remeber easily at the columns of the matrix. 

    % NumMec = 1;
    % NumSeq = 2;
    % NumImg = 3;
    % stampN = 4;
    % numName = 5;
    % numAdd = 6;
    % door_number = 7
    % D1 = 8;
    % D2 = 9;
    % D3 = 10;
    % D4 = 11;
    % D5 = 12;
    % D6 = 13;
    % D7 = 14;
    % D8 = 15;


    statistics_matrix = zeros(rows, columns);

    %I already know that in the first column i have to put my student
    %number so i'm gonna do it right now
    %same thing for NumSeq and NumImg
    
    for i = 1:size(images_in_previous_directory)

        %adding the NumMec
        statistics_matrix(i, 1) = 122422;

        %finding the num_seq and num_img
        num_seq_num_img = '(\d+)_+(\d+)\.png';
        found = regexp(images_in_previous_directory(i).name, num_seq_num_img, "tokens");
        %adding NumSeq & NumImg

        if ~isempty(found)
            statistics_matrix(i, 2) = str2double(found{1}{1});
            statistics_matrix(i, 3) = str2double(found{1}{2});
        end
    end

    %Now for all the images we call a function which will handle all the
    %images

    for i = 1:size(list_of_all_images, 2)
       statistics_matrix = start_processing(list_of_all_images{i}, statistics_matrix, i);
    end

    

    %creating the txt file

    file_of_results = fopen('tp2_122422.txt', 'w');

    for i = 1:size(statistics_matrix, 1)

        result_row = ' ';

        for j = 1:size(statistics_matrix, 2)

            result_row = [result_row, num2str(statistics_matrix(i, j))];

            if j < size(statistics_matrix, 2)

                result_row = [result_row, ',']; 

            end

        end
        fprintf(file_of_results, '%s\n', result_row);
    end

    fclose(file_of_results);
    % 
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

function output_image = filtered_image(input_image)

    % this operation is useful to understand if we need to do a more
    % complex operation on the image while filtering or smìimpky binarize
    % it
    
    gray_image = rgb2gray(input_image);
    adjusted_image = imadjust(gray_image);
    output_image = ~imbinarize(adjusted_image);
end


function new_statistics_matrix = start_processing(image, statistics_matrix, index)

    %In this function we will start the process of the images. This will be
    %the scheme:

    % 1. Filtering the image
  
    images_filtered = filtered_image(image);

    if statistics_matrix(index, 2) == 120 %i rotate the image only if is a 120 image
        images_filtered = rotate(images_filtered);
    end

    num_seq = statistics_matrix(index, 2);
    
    %find number of word in the name and number of word in the address door
    %number and postal code
    [num_of_word_in_name, num_of_word_in_address, door, pc] = count_word_in_name_and_address(images_filtered, num_seq, index);

    statistics_matrix(index, 5) = statistics_matrix(index, 5) + num_of_word_in_name;
    statistics_matrix(index, 6) = statistics_matrix(index, 6) + num_of_word_in_address; 
    statistics_matrix(index, 7) = statistics_matrix(index, 7) + str2num(door);
    
    for j = 1:numel(pc)
        statistics_matrix(index, 7+j) = statistics_matrix(index, 7+j) + (int32(pc{j})-1);
    end

    %we recognize the selo

    [type_of_selo] = classify_selos(image, num_seq, index);

    statistics_matrix(index, 4) = statistics_matrix(index, 4) + type_of_selo;

    %we update the matrix
    new_statistics_matrix = statistics_matrix;

end

function rotated_image = rotate(image)

    edges = edge(image, 'canny');
    [H, theta, rho] = hough(edges);

    % finding the peaks of the hough transform
    peaks = houghpeaks(H, 5, 'threshold', ceil(0.3*max(H(:))));
    
    % Extracting lines
    lines = houghlines(edges, theta, rho, peaks, 'FillGap', 20, 'MinLength', 40);
        
    angles = [lines.theta];
    if isempty(angles)
        rotationAngle = 0;  
    else
        rotationAngle = median(angles);
    end
    
    % computing the angle
    if rotationAngle > 45
        rotationAngle = rotationAngle - 90;
    elseif rotationAngle < -45
        rotationAngle = rotationAngle + 90;
    end

    % rotating the image
    rotated_image = imrotate(image, rotationAngle);
    
end

function cropped_image = cut_image_100(image)
    
    [rows, cols, ~] = size(image);
    margin = 0.25; 
    x_start = round(cols * margin);
    y_start = round(rows * margin);
    x_end = cols - round(cols * margin);
    y_end = rows - round(rows * margin);

    
    cropped_image = image(y_start:y_end+30, x_start-30:x_end - 30, :);

end

function cropped_image = cut_image_120(image)
    
    [rows, cols, ~] = size(image);
    margin = 0.25; 
    x_start = round(cols * margin);
    y_start = round(rows * margin);
    x_end = cols - round(cols * margin);
    y_end = rows - round(rows * margin);

    
    cropped_image = image(y_start + 60 :y_end+30, x_start-50:x_end - 90, :);

end

function [number_name, number_address, door, pc ] = count_word_in_name_and_address(image, num_seq, index)
    
%This function returns the number of words that are in each name and each
%address

    % First of all we delete the borders of the image
    if(num_seq == 100)
        cropped_image_non_final = cut_image_100(image);
         cropped_image_non_final = imclearborder( cropped_image_non_final);
    else
        cropped_image_non_final= cut_image_120(image);
         cropped_image_non_final = imclearborder( cropped_image_non_final);
    end
    
    cropped_image = bwareaopen( cropped_image_non_final, 45);
    
    %now we detect the text
    words_in_postmail = ocr(cropped_image);

    textLineBoundingBoxes = words_in_postmail.TextLineBoundingBoxes;
    wordBoundingBoxes = words_in_postmail.WordBoundingBoxes;
    words_in_postmail.TextLines{1};

    initial_row = 1;
    name = words_in_postmail.TextLines{initial_row };
    address = words_in_postmail.TextLines{initial_row  + 1};

    %we extract the token and then count
    token_name = strsplit(strtrim(name));
    token_address = strsplit(strtrim(address));

    number_name = length(token_name);
    number_address = length(token_address) - 1;

    %i extract the postal code and the door number

    thirdLineBox = textLineBoundingBoxes(initial_row +2, :);
    thirdLineBox(1) = thirdLineBox(1) - 30;
    thirdLineBox(3) = thirdLineBox(3) + 70;
    second_row = textLineBoundingBoxes(initial_row + 1, :);
    secondline = imcrop(cropped_image_non_final, second_row);

    
    %the only way to extract the number correctly seems to be using 
    %structural element

    structural_element = strel('square', 10);
    second_line = imclose(secondline, structural_element);
    
    [label_seocond_line, number_of_label] = bwlabel(second_line);

    properties_second_line = regionprops(label_seocond_line, 'BoundingBox');

    boundin_boxes = cat(1, properties_second_line.BoundingBox);

    max_x = - inf;
    last_bounding_box = [];
    for i = 1:number_of_label
        current_bounding_box = boundin_boxes(i, :);
        if current_bounding_box(1) > max_x
            max_x = current_bounding_box(1);
            last_bounding_box = current_bounding_box;
        end

    end
   
    %extracting the number for the recognition
    postal_code = imcrop(cropped_image_non_final, thirdLineBox);
    door_number = imcrop(secondline, last_bounding_box);
  
    [door_number_recognition, postal_code_recognition] = recognize_numbers(door_number, postal_code, index);
    door = door_number_recognition;
    pc = postal_code_recognition;

end

function [door_number_recognition, postal_code_recognition] = recognize_numbers(door_number, postal_code, index)

    load("mnist_svpi2024_tp2_122422.mat", 'net');
    door_number_recognition = '';
    postal_code_recognition = {};

    numbers_door = {}; 
    numbers_pc = {};

    door_number = bwareaopen(door_number, 5);
    [door_label, door_num] = bwlabel(door_number);
    door_bounding_box = regionprops(door_label, 'BoundingBox');

    for i = 1: door_num
            current_number = imcrop(door_number, door_bounding_box(i).BoundingBox);
            numbers_door{end+1} = current_number;
    end

    structured_element = strel('disk', 5);
    postal_code_structural_element = imclose(postal_code, structured_element);
    [pc_label, pc_num] = bwlabel(postal_code_structural_element);
    pc_bounding_box = regionprops(pc_label, 'BoundingBox');

    for i = 1: pc_num
        if pc_bounding_box(i).BoundingBox(3) * pc_bounding_box(i).BoundingBox(4) > 10
            current_number = imcrop(postal_code, pc_bounding_box(i).BoundingBox);
            numbers_pc{end+1} = current_number;
        end
    end

    %to better clssify the number i resize them to have a 28x28 pixel
    %dimension
    for i = 1:numel(numbers_door)
        current_number = numbers_door{i};
        [row, column] = size(current_number);
        if row > 28 || column > 28
            current_number = imresize(current_number, [28 28]);
            row = 28;
            column = 28;
        end
        padded_image = zeros(28, 28);
        row_offset = floor((28 - row) / 2) + 1;
        col_offset = floor((28 - column) / 2) + 1;
        end_row_idx = min(row_offset + row - 1, 28);
        end_col_idx = min(col_offset + column - 1, 28);
        padded_image(row_offset:end_row_idx, col_offset:end_col_idx) = current_number;
        result = classify(net, padded_image);
        door_number_recognition = strcat(door_number_recognition, char(result));
    end

    %postal code recognition
    for i = 1:numel(numbers_pc)
        current_number = numbers_pc{i};
        [row, column] = size(current_number);
        if row > 28 || column > 28
            current_number = imresize(current_number, [28 28]);
            row = 28;
            column = 28;
        end
        padded_image = zeros(28, 28);
        row_offset = floor((28 - row) / 2) + 1;
        col_offset = floor((28 - column) / 2) + 1;
        end_row_idx = min(row_offset + row - 1, 28);
        end_col_idx = min(col_offset + column - 1, 28);
        padded_image(row_offset:end_row_idx, col_offset:end_col_idx) = current_number;
        result = classify(net, padded_image);
        postal_code_recognition{end+1} = result;
    end
end

function [type_of_selo] = classify_selos(image, num_seq, index)

    %This function take the image and try to classify the selo that is on
    %the image

    type_of_selo = 0; %at the moment we assume there is no selo

                        %ffa       sol       ecc    ca(convex area)
    selos_properties = [0.4554    0.9548    0.7046  189.5300;        %type 1
                        0.4825    0.9575    0.5883  215.4500;        %type 2
                        0.4770    0.9587    0.7451  200.6700;        %type 3
                        0.4614    0.9545    0.6272  172.8300;        %type 4
                        0.4577    0.9479    0.4825  154.2200;        %type 5
                        0.4656    0.9584    0.5871  215.4600;        %type 6
                        0.4724    0.9601    0.5114  228.4100];       %type 7

    p1 = selos_properties(1, :);
    p2 = selos_properties(2, :);
    p3 = selos_properties(3, :);
    p4 = selos_properties(4, :);
    p5 = selos_properties(5, :);
    p6 = selos_properties(6, :);
    p7 = selos_properties(7, :);

    %The light blue postmail are giving me problem so i change the bg in
    %black
    backgroundMask = image(:,:,1) == 194 & ...
                 image(:,:,2) == 227 & ...
                 image(:,:,3) == 229;
    image(repmat(backgroundMask, [1 1 3])) = 0;
    [rows, ~] = size(image);

    if (num_seq == 120)
        selos_pos_image= image(50:(rows/2)+50, end-360:end-60, :);
    else
        selos_pos_image= image(60:(rows/2)+40, end-315:end-60, :);
    end

    %doing some filtering to a better classification
    selos_pos_gray = rgb2gray(selos_pos_image);
    selos_pos_adj = imadjust(selos_pos_gray);
    white_mask = (selos_pos_adj  == 255 |  selos_pos_adj  > 240 | all(selos_pos_adj >= 0.9, 2));
    selos_pos_adj(white_mask) = 0;
    selos_pos_bin = imbinarize(selos_pos_adj);
    selos_pos_bin = imclearborder(selos_pos_bin);
    selos_pos_bin = imfill(selos_pos_bin, 'holes');
    selos_pos_bin = imclearborder(selos_pos_bin);
    selos_pos_bin = bwareaopen(selos_pos_bin, 2);
    se = strel('diamond', 4);
    selos_pos_se = imclose(selos_pos_bin, se);
    selos_pos_filled = imfill(selos_pos_se, 'holes');
    selos_pos_filled = bwareaopen(selos_pos_filled, 1200);
    selos_pos_filled = imclearborder(selos_pos_filled); 
    [L_selos, N] = bwlabel(selos_pos_filled);

    if N == 0 %no regions found
        type_of_selo = 0;
        return
    end

    %computing current selo properties
    selos_pos_props = regionprops(L_selos, 'Solidity', 'Eccentricity', 'Circularity', 'ConvexArea'); 
    sol_selos = [selos_pos_props.Solidity]';
    ecc_selos = [selos_pos_props.Eccentricity]';
    ffa_selos = [selos_pos_props.Circularity]';
    ca_selos = [selos_pos_props.ConvexArea]'/100;
    Patts_selos = [ffa_selos sol_selos ecc_selos ca_selos];

    %computing the distances
    d1 = vecnorm(Patts_selos-p1, 2, 2); 
    d2 = vecnorm(Patts_selos-p2, 2, 2);
    d3 = vecnorm(Patts_selos-p3, 2, 2);
    d4 = vecnorm(Patts_selos-p4, 2, 2);
    d5 = vecnorm(Patts_selos-p5, 2, 2);
    d6 = vecnorm(Patts_selos-p6, 2, 2);
    d7 = vecnorm(Patts_selos-p7, 2, 2);

    distances = [d1 d2 d3 d4 d5 d6 d7];
    
    [~, i] = min(distances);

    %it shold be only one region detected
    if(numel(i) > 1)
        type_of_selo = 0;
        return
    end

    % Sometimes type 6 and 2 are confused. So if in the selo there is
    % yellow it means that is the type = 2
    if i == 6 && ~any((selos_pos_image(:,:,1) <= 161 & selos_pos_image(:,:,1) >= 158) & ...
                 selos_pos_image(:,:,2) == 91 & ...
                 (selos_pos_image(:,:,3) >= 84 & selos_pos_image(:,:,3) <= 88), 'all')
        i = 2;
    end

    type_of_selo = i;
end



