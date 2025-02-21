from PIL import Image
import os
import io
from natsort import natsorted
from reportlab.lib.pagesizes import letter, landscape
from reportlab.pdfgen import canvas

# path information
project_path = "2024-06-13-tp53_competition_sameMut_sameP"
path = f"/net/gs/vol1/home/huntc10/proj/HomeostaticEpithelium/{project_path}/results/"

## take 10 frames from directory and stitch them together into 1 sequence image
def combine_images_in_directory(directory_path, output_path):
    # Get list of files in directory
    files = os.listdir(directory_path)
    prefixes = ["frame_1_", "frame_2_", "frame_3_", "frame_4_", "frame_5_", "frame_6_", "frame_7_", "frame_8_", "frame_9_", "frame_10_"]
    images = [file for file in files if file.endswith(".png") and any(file.startswith(prefix) for prefix in prefixes)]

    images = natsorted(images) # Sort images by filename
    
    # Open images and get their dimensions
    image_objects = [Image.open(os.path.join(directory_path, img)) for img in images]
    widths, heights = zip(*(img.size for img in image_objects))

    # Calculate total width and height for the final image
    total_width = sum(widths)
    max_height = max(heights)

    # Create a new blank image with the calculated dimensions
    combined_image = Image.new("RGB", (total_width, max_height))

    # Paste each image into the combined image
    x_offset = 0
    for img in image_objects:
        combined_image.paste(img, (x_offset, 0))
        x_offset += img.width

    combined_image.save(output_path) # Save the combined image

## use stitched images to make a pdf file of all image sequences
def create_pdf(output_directory, directory_paths, block):
    c = canvas.Canvas(f'{output_directory}summary_timeseries_block_{block}.pdf', pagesize=landscape(letter))
    y_offset = 0
    x_offset = 0
    for directory in directory_paths:
        files = os.listdir(directory)
        combined_image = [file for file in files if file.startswith("all_combined_block_")]

        if(len(combined_image) >0):
            #print(combined_image[0])
            combined_image_path = f'{directory}{combined_image[0]}'
            img = Image.open(combined_image_path)
            width, height = img.size
            aspect_ratio = width / height
            max_width = 790
            max_height = 79
            # Determine scaling factor to fit within max_width and max_height
            width_scale = max_width / width
            height_scale = max_height / height
            scale_factor = min(width_scale, height_scale, 1.0)  # Ensure the image is scaled down, not up
            # Scale the image
            new_width = int(width * scale_factor)
            new_height = int(height * scale_factor)
            img = img.resize((new_width, new_height), Image.LANCZOS)

            if y_offset + new_height > landscape(letter)[1]:
            # Add a new page
                c.showPage()
                y_offset = 0

            # Draw the image on the PDF
            c.drawInlineImage(img, 0, y_offset, width=new_width, height=new_height)

            # Update the y-coordinate for the next image
            y_offset += new_height
            y_offset += 5

            # Add some space between images
            x_offset += 0

    c.save()


############
###### Main
############
    
# blocking prob and replicate information
corrBlocks = [0,0.01,0.1]
replicates = range(0,100)
# generate filepaths for each directory with the set of images

for k in corrBlocks:
    directory_paths = []
    for j in replicates:
        frame_folder = f'{path}/animation_frames/block_{k}_rep_{j}/'
        if (os.path.exists(frame_folder)):
            directory_paths.append(frame_folder)

        # loop through each directory of images and combine the 10 images
    for directory in directory_paths:
        output_filename = f"{directory}all_combined_block_{k}.png"
        # Combine images in directory and save to output directory
        combine_images_in_directory(directory, output_filename)
        print(f"all images combined for block: {k} in {directory}")
    print(f"creating pdf for block: {k}")
    pdf_directory = f"{path}"
    create_pdf(pdf_directory, directory_paths, k)
    print(f"created pdf for block: {k}")