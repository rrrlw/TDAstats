library("hexSticker")

imgurl <- "barcode.png"
sticker(imgurl, package = "TDAstats",
        p_size = 20, p_y = 1.45,
        s_x = 1, s_y = 0.75, s_width = 0.6,
        h_fill = "black", h_color = "white",
        filename = "TDAstatsSticker.png")
