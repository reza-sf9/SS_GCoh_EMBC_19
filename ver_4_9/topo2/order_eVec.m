function [order_eVec] = order_eVec(eVec, order_num)
% Read data_layout_ds24_easycapm25.docx

% input data is organized based dataset (DSI-24)
% I should change the order based the layout (easycapM25)

switch order_num
    
    % for Ali's data set
    case 64
        order_dataset = [18, 17, 48, 54, 24, 12, 11, 41, 19, 13, 6, 5, 2, 36, 37, 42, 49, 60, 55, 43, 38, 33, 1, 3, 7, 20, 25, 14,...
            8, 4, 9, 35, 34, 39, 44, 50, 61, 56, 51, 45, 40, 10, 16, 15, 21, 26, 29, 30, 27, 22, 23, 47, 53, 46, 52, 57, 62,...
            58, 59, 28, 31, 32, 64, 63];
        
    % for Ali's data set
    case 32
        order_dataset = [12, 11, 24, 28, 13, 9, 8, 4, 3, 1, 18, 25, 31, 29, 2, 5, 6, ...
            17, 19, 21, 26, 32, 30, 22, 20, 7, 10, 16, 14, 23, 27, 15];
        
        
        % for Yalda's data set 
    case 20
        order_dataset = [10, 16, 11, 17, 3, 4, 5, 20, 6, 8, 2, 12, 18, 19, 7, 9, 1, 13, 14, 15];
end

m = size(eVec);
order_eVec = zeros(m(1), m(2));

for i=1 : m(1)
    order_eVec(i, 1) = eVec(order_dataset(1, i), 1);
end

end