clc; clear variables; close all; 

lot_number = randperm(20, 6); 
N = 10e6; 
matches = 0; 
for n = 1:N
    ticket = randperm(20, 6); 
    if num_match(ticket, lot_number) == 6
        matches = matches + 1; 
    end

end

win_prob = matches / N; 
