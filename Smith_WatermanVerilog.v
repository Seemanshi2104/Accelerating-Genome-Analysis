`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 23.03.2024 15:42:00
// Design Name: 
// Module Name: Smith_Waterman
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////

module Smith_Waterman(    
    input [0:15] sequence1,
    input [0:13] sequence2,
    output reg [7:0] Ldistance);

parameter SEQ1_LENGTH = 16;
parameter SEQ2_LENGTH = 14;
parameter MATCH_REWARD = 2;
parameter MISMATCH_PENALTY = 1;
parameter GAP_PENALTY = 0;

reg [3:0] main [0:SEQ1_LENGTH/2][0:SEQ2_LENGTH/2];

integer i, j;
// Initialize first row and first column of main as 0
initial begin
    main[0][0] = 0;
    for (i = 1; i <= SEQ1_LENGTH/2; i = i + 1)
        main[i][0] = 0;
    for (j = 1; j <= SEQ2_LENGTH/2; j = j + 1)
        main[0][j] = 0;
end

reg [3:0] value_from_left, value_from_right, value_from_diag, check;
always @(*) begin
    for (i = 1; i <= SEQ1_LENGTH/2; i = i + 1) begin
        for (j = 1; j <= SEQ2_LENGTH/2; j = j + 1) begin

            value_from_left = main[i][j-1] + GAP_PENALTY;
            value_from_right = main[i-1][j] + GAP_PENALTY;
            if (sequence1[2*i-1] == sequence2[2*j-1] && sequence1[2*i-2] == sequence2[2*j-2])
                check = MATCH_REWARD;
            else
                check = MISMATCH_PENALTY;
            value_from_diag = main[i-1][j-1] + check;
            main[i][j] = value_from_diag > value_from_left ? (value_from_diag > value_from_right ? value_from_diag : value_from_right) : (value_from_left > value_from_right ? value_from_left : value_from_right);
        end
    end
end

reg [3:0] max_value;
integer x, y;
integer p, q;

always @(*) begin
    max_value = main[1][1]; 
    for (p = 1; p <= SEQ1_LENGTH/2; p = p + 1) begin
        for (q = 1; q <= SEQ2_LENGTH/2; q = q + 1) begin
            if (main[p][q] > max_value) begin
                max_value = main[p][q]; 
                x = p; 
                y = q;
            end
        end
    end
end

//reg [3:0] match_checker_matrix[0:(SEQ1_LENGTH/2)-1][0:(SEQ2_LENGTH/2)-1];

//always @(*)begin

//    for (i = 0; i <= SEQ1_LENGTH/2; i = i + 1) begin
//        for (j = 0; j <= SEQ2_LENGTH/2; j = j + 1) begin
//            if (sequence1[2*i] == sequence2[2*j] && sequence1[2*i+1] == sequence2[2*j+1])
//                match_checker_matrix[i][j] = MATCH_REWARD;
//            else
//                match_checker_matrix[i][j] = MISMATCH_PENALTY;
//        end
//    end

//end

reg [0:SEQ1_LENGTH+1]aligned1;
reg [0:SEQ1_LENGTH+1]aligned2;
reg flag;
initial flag=0;
reg [3:0] temp;

initial Ldistance = 0; // Ldistance represents the Levensthtein distance, which counts the number of mismatch and gaps

always @(*)begin

    while (main[x][y] != 0) begin
        temp = (sequence1[2*x-2] == sequence2[2*y-2]) && (sequence1[2*x-1] == sequence2[2*y-1]);
        if (main[x][y] == main[x-1][y-1] + (temp ? MATCH_REWARD : MISMATCH_PENALTY)) begin
                    aligned1[1] = sequence1[2*x-1]; 
                    aligned1[0] = sequence1[2*x-2];
                    
                    aligned2[1] = sequence2[2*y-1]; 
                    aligned2[0] = sequence2[2*y-2];
                    
                    if ((x != 1) || (y != 1))begin
                        aligned1 = aligned1 >> 2;
                        aligned2 = aligned2 >> 2;                
                     end
                     
                     if (temp == 0)
                        Ldistance = Ldistance + 1;
            x = x - 1;
            y = y - 1;

        end
        else if (main[x][y] == (main[x-1][y] + GAP_PENALTY)) begin
            aligned2[1] = 1'bX;
            aligned2[0] = 1'bX;
            
            aligned1[1] = sequence1[2*x-1]; 
            aligned1[0] = sequence1[2*x-2];
            
            if ((x != 1) || (y != 1))begin
                aligned1 = aligned1 >> 2;
                aligned2 = aligned2 >> 2;                
            end
             
            Ldistance = Ldistance + 1;
            x = x - 1;
        end
        else begin
            aligned1[1] = 1'bX;
            aligned1[0] = 1'bX;
            
            aligned2[1] = sequence2[2*y-1]; 
            aligned2[0] = sequence2[2*y-2];
            
            if ((x != 1) || (y != 1))begin
                aligned1 = aligned1 >> 2;
                aligned2 = aligned2 >> 2;                
            end
            
            Ldistance = Ldistance + 1;
            y = y - 1;
        end
    end
end
endmodule
