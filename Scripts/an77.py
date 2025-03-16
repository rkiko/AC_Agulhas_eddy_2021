import fitz
from pathlib import Path
home = str(Path.home())

file1="%s/GIT/AC_Agulhas_eddy_2021/Plots/Fig_Main_v08/Fig03a_v08.pdf" % (home)
file2="%s/GIT/AC_Agulhas_eddy_2021/Plots/Fig_Main_v08/Fig03b_v08.pdf" % (home)

# Open the PDF files
pdf1 = fitz.open(file1)  # The PDF to overlay (Fig03a)
pdf2 = fitz.open(file2)  # The base PDF (Fig03b)

# Get the first page of each PDF
page1 = pdf1.load_page(0)  # First page of Fig03a
page2 = pdf2.load_page(0)  # First page of Fig03b

#page dimension
width1, height1 = page1.rect.width, page1.rect.height  # Dimensions of Fig03a
width2, height2 = page2.rect.width, page2.rect.height  # Dimensions of Fig03b

# Scale Fig03a to match the width of Fig03b while maintaining aspect ratio
scale_factor = width2 / width1  # Scaling factor
scaled_height1 = height1 * scale_factor  # New height of Fig03a after scaling

# Calculate the total height of the new PDF
total_height = scaled_height1 + height2
max_width = width2  # Use the width of Fig03b as the page width

# Create a new PDF to store the result
output_pdf = fitz.open()

# Create a new page with the combined height and width of Fig03b
new_page = output_pdf.new_page(width=max_width, height=total_height)

# Place Fig03a at the top (scaled to match the width of Fig03b)
a=0.8
new_page.show_pdf_page(
    fitz.Rect(30, 30, width2*a, scaled_height1*a+30),  # Position at the top
    pdf1,  # PDF to insert (Fig03a)
    0,     # Page number
)

# Place Fig03b at the bottom
new_page.show_pdf_page(
    fitz.Rect(0, scaled_height1, width2, scaled_height1 + height2),  # Position at the bottom
    pdf2,  # PDF to insert (Fig03b)
    0,     # Page number
)

# # Define the position and size of the overlay (Fig03a)
# # Adjust these values to change the position and size
# x0 = 50  # Left margin
# y0 = 150  # Bottom margin
# x1 = page2.rect.width - 50  # Right margin
# y1 = page2.rect.height - 10  # Top margin
#
# overlay_rect = fitz.Rect(x0, y0, x1, y1)
#
# # Overlay the second PDF (Fig03a) onto the new page
# new_page.show_pdf_page(overlay_rect, pdf1, 0)

# Save the result
output_pdf.save("%s/GIT/AC_Agulhas_eddy_2021/Plots/Fig_Main_v08/merged_output.pdf" % (home))
output_pdf.close()