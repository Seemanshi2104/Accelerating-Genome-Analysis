module Murmurhash (
    input wire [41:0] KEY,
    output reg [31:0] EI
);

    reg [41:0] D;

    always @* begin
        D = KEY ^ (KEY >> 33);
    end

    always @* begin
        D = D * 32'hFF51AFD7ED558CCD;
    end

    always @* begin
        D = D ^ (KEY >> 33);
    end

    always @* begin
        D = D * 32'hC4CEB9FE1A85EC53;
    end

    always @* begin
        D = D ^ (KEY >> 33);
    end

    always @* begin
        EI = D[31:0];
    end

endmodule
