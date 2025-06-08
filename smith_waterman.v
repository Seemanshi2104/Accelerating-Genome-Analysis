module smith_waterman #(
    parameter REF_LEN = 15,
    parameter QUERY_LEN = 10,
    parameter BASE_WIDTH = 2,
    parameter ALIGN_LEN = REF_LEN + QUERY_LEN // max possible alignment length
)(
    input clk,
    input rst,
    input [REF_LEN*BASE_WIDTH-1:0] ref_seq,
    input [QUERY_LEN*BASE_WIDTH-1:0] query_seq,
    output reg [ALIGN_LEN*BASE_WIDTH-1:0] aligned_ref_seq,
    output reg [ALIGN_LEN*BASE_WIDTH-1:0] aligned_query_seq,
    output reg [7:0] alignment_length  // number of valid aligned pairs
);

    localparam MATCH_SCORE    =  3;
    localparam MISMATCH_SCORE = -1;
    localparam GAP_SCORE      = -2;

    integer H[0:REF_LEN][0:QUERY_LEN];
    integer traceback[0:REF_LEN][0:QUERY_LEN];

    integer i, j, ti, tj;
    integer score_diag, score_up, score_left, score;
    integer max_score, max_i, max_j;
    reg [1:0] ref_base, query_base;

    // Helpers to extract bases
    function [1:0] get_base(input [REF_LEN*BASE_WIDTH-1:0] seq, input integer index);
        get_base = seq[(REF_LEN-1-index)*BASE_WIDTH +: BASE_WIDTH];
    endfunction

    function [1:0] get_query_base(input [QUERY_LEN*BASE_WIDTH-1:0] seq, input integer index);
        get_query_base = seq[(QUERY_LEN-1-index)*BASE_WIDTH +: BASE_WIDTH];
    endfunction

    // Main logic
    always @(posedge clk or posedge rst) begin
        if (rst) begin
            // Reset everything
            for (i = 0; i <= REF_LEN; i = i + 1)
                for (j = 0; j <= QUERY_LEN; j = j + 1) begin
                    H[i][j] = 0;
                    traceback[i][j] = 0;
                end
            max_score = 0;
            max_i = 0;
            max_j = 0;
            aligned_ref_seq = 0;
            aligned_query_seq = 0;
            alignment_length = 0;
        end else begin
            // Fill score and traceback matrices
            for (i = 1; i <= REF_LEN; i = i + 1) begin
                for (j = 1; j <= QUERY_LEN; j = j + 1) begin
                    ref_base = get_base(ref_seq, i-1);
                    query_base = get_query_base(query_seq, j-1);

                    if (ref_base == query_base)
                        score_diag = H[i-1][j-1] + MATCH_SCORE;
                    else
                        score_diag = H[i-1][j-1] + MISMATCH_SCORE;

                    score_up   = H[i-1][j] + GAP_SCORE;
                    score_left = H[i][j-1] + GAP_SCORE;

                    score = score_diag;
                    traceback[i][j] = 1; // diagonal

                    if (score_up > score) begin
                        score = score_up;
                        traceback[i][j] = 2; // up
                    end
                    if (score_left > score) begin
                        score = score_left;
                        traceback[i][j] = 3; // left
                    end
                    if (score < 0) begin
                        score = 0;
                        traceback[i][j] = 0;
                    end

                    H[i][j] = score;

                    if (score > max_score) begin
                        max_score = score;
                        max_i = i;
                        max_j = j;
                    end
                end
            end

            // Traceback
            ti = max_i;
            tj = max_j;
            alignment_length = 0;
            aligned_ref_seq = 0;
            aligned_query_seq = 0;

            while (ti > 0 && tj > 0 && H[ti][tj] > 0 && alignment_length < ALIGN_LEN) begin
                case (traceback[ti][tj])
                    1: begin // diagonal (match/mismatch)
                        aligned_ref_seq[alignment_length*BASE_WIDTH +: BASE_WIDTH] = get_base(ref_seq, ti-1);
                        aligned_query_seq[alignment_length*BASE_WIDTH +: BASE_WIDTH] = get_query_base(query_seq, tj-1);
                        ti = ti - 1;
                        tj = tj - 1;
                    end
                    2: begin // up (gap in query)
                        aligned_ref_seq[alignment_length*BASE_WIDTH +: BASE_WIDTH] = get_base(ref_seq, ti-1);
                        aligned_query_seq[alignment_length*BASE_WIDTH +: BASE_WIDTH] = 2'bxx; // GAP
                        ti = ti - 1;
                    end
                    3: begin // left (gap in ref)
                        aligned_ref_seq[alignment_length*BASE_WIDTH +: BASE_WIDTH] = 2'bxx; // GAP
                        aligned_query_seq[alignment_length*BASE_WIDTH +: BASE_WIDTH] = get_query_base(query_seq, tj-1);
                        tj = tj - 1;
                    end
                    default: begin
                        ti = 0;
                        tj = 0;
                    end
                endcase
                alignment_length = alignment_length + 1;
            end
        end
    end
endmodule








module tb_smith_waterman;

    // Parameters
    localparam REF_LEN = 15;
    localparam QUERY_LEN = 10;
    localparam BASE_WIDTH = 2;
    localparam ALIGN_LEN = REF_LEN + QUERY_LEN;

    // Signals
    reg clk;
    reg rst;
    reg [REF_LEN*BASE_WIDTH-1:0] ref_seq;
    reg [QUERY_LEN*BASE_WIDTH-1:0] query_seq;
    wire [ALIGN_LEN*BASE_WIDTH-1:0] aligned_ref_seq;
    wire [ALIGN_LEN*BASE_WIDTH-1:0] aligned_query_seq;
    wire [7:0] alignment_length;

    // Instantiate the Smith-Waterman module
    smith_waterman #(
        .REF_LEN(REF_LEN),
        .QUERY_LEN(QUERY_LEN),
        .BASE_WIDTH(BASE_WIDTH),
        .ALIGN_LEN(ALIGN_LEN)
    ) uut (
        .clk(clk),
        .rst(rst),
        .ref_seq(ref_seq),
        .query_seq(query_seq),
        .aligned_ref_seq(aligned_ref_seq),
        .aligned_query_seq(aligned_query_seq),
        .alignment_length(alignment_length)
    );

    // Clock generation
    always begin
        #5 clk = ~clk;  // 10ns period
    end

    // Function to decode 2-bit base to character
    function [7:0] base_to_char(input [1:0] base);
        case (base)
            2'b00: base_to_char = "A";
            2'b01: base_to_char = "T";
            2'b10: base_to_char = "G";
            2'b11: base_to_char = "C"; 
            2'bxx: base_to_char = "-"; 
            default: base_to_char = "?";
        endcase
    endfunction

    // Task to print aligned sequences
    task print_alignment;
        integer i;
        reg [1:0] ref_b, query_b;
        begin
            $write("Aligned Reference: ");
            for (i = alignment_length-1; i >= 0; i = i - 1) begin
                ref_b = aligned_ref_seq[i*2 +: 2];
                $write("%s", base_to_char(ref_b));
            end
            $write("\n");

            $write("Aligned Query   : ");
            for (i = alignment_length-1; i >= 0; i = i - 1) begin
                query_b = aligned_query_seq[i*2 +: 2];
                $write("%s", base_to_char(query_b));
            end
            $write("\n");
        end
    endtask

    // Stimulus
    initial begin
        // Initialize
        clk = 0;
        rst = 1;
        ref_seq = 0;
        query_seq = 0;

        // Apply reset
        #10;
        rst = 0;

        // Reference and query (ACGTACGTACGTA and ACGTACGTA)
        ref_seq = {2'b10 , 2'b01, 2'b00, 2'b01, 2'b10, 2'b11, 2'b00, 2'b01, 2'b01, 2'b10, 2'b11, 2'b00, 2'b01, 2'b10, 2'b10};
        query_seq = {2'b00, 2'b01, 2'b10, 2'b11, 2'b00, 2'b01, 2'b01, 2'b10,  2'b11,2'b00};

        // Wait for alignment to compute
        #200;

        // Show result
        $display("\n>> Smith-Waterman Local Alignment Result:");
        print_alignment();

        // End simulation
        $finish;
    end

endmodule





