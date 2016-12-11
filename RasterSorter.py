"""
Sorts raster images in directory by features in given shapefile
"""

import gdal
import os
from osgeo import ogr
import shapefile
import shutil
import ntpath
import datetime
from hashlib import sha1
import logging

class Raster:
    def __init__(self):
        pass

    def open_raster(self, raster_file):
        """
        Opens the raster data using gdal
        :param raster_file: input file
        :return:
        """

        geotransform = None
        data = None
        columns = None
        rows = None

        logging.info("Opening raster file {0}".format(raster_file))

        try:
            data = gdal.Open(raster_file)
        except Exception as e:
            logging.error("\n***** Could not load raster file {0]. The tool return {1}".format(raster_file, e))


        logging.info("Getting geotransform of raster file {0}".format(raster_file))

        try:
            geotransform = data.GetGeoTransform()
        except Exception as e:
            logging.error("\n***** Could not get geotransform of raster file {0}. The tool returned {1}"
                          .format(raster_file, e))

        if data:
            logging.info("Getting rows and columns of raster file {0}".format(raster_file))
            columns = data.RasterXSize
            rows = data.RasterYSize

        return geotransform, columns, rows

    def get_raster_extents(self, raster_file):
        """
        Gets raster extents
        :param raster_file: input raster
        :return:
        """

        ext = []
        geotransform, columns, rows = self.open_raster(raster_file)

        if geotransform and columns and rows:
            x_array = [0, columns]
            y_array = [0, rows]

            for px in x_array:
                for py in y_array:
                    x = geotransform[0] + (px * geotransform[1]) + (py * geotransform[2])
                    y = geotransform[3] + (px * geotransform[4]) + (py * geotransform[5])
                y_array.reverse()
        return ext

    def generate_raster_bbox(self, raster_file):
        """
        Generate bounding box for raster
        :param raster_file: raster input file
        :return:
        """

        raster = gdal.Open(raster_file)
        transform = raster.GetGeoTransform()
        pixel_width = transform[1]
        pixel_height = transform[5]
        cols = raster.RasterXSize
        rows = raster.RasterYSize

        xleft = transform[0]
        ytop = transform[3]
        xright = xleft + (cols * pixel_width)
        ybottom = ytop - (rows * pixel_height)

        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(xleft, ytop)
        ring.AddPoint(xleft, ybottom)
        ring.AddPoint(xright, ytop)
        ring.AddPoint(xright, ybottom)
        ring.AddPoint(xleft, ytop)

        raster_geometry = ogr.Geometry(ogr.wkbPolygon)
        raster_geometry.AddGeometry(ring)

        return raster_geometry


class Shapes:
    def __init__(self):
        pass

    def get_vector_geometry(self, vector_data):
        """
        Generates vector geometry
        :param vector_data: vector data input
        :return:
        """

        geometries = []
        vector = gdal.Open(vector_data)
        layer = vector.getLayer()
        n_features = vector.GetFeatureCount()

        for feature_number in range(n_features):
            feature = layer.GetFeature(feature_number)
            vector_geometry = feature.GetGeometryRef()
            geometries.append(vector_geometry)

        return geometries

    def get_feature_bounds(self, shape_data):
        """
        Gets the bounding coordinates of each feature in shapefile
        :param shape_data: Shapefile input
        :return:
        """

        shp = shapefile.Reader(shape_data)
        shapes = shp.shapes()
        records = shp.records()
        n_shapes = len(list(shp.iterShapes()))

        results = []

        for shp_counter in range(n_shapes):
            # get bounding box
            bbox = shapes[shp_counter].bbox
            subsection_name = records[shp_counter[9]]  #TODO: Get user input for column to sort images on

            results.append((subsection_name, bbox))

        return results


def make_directory(directory):
    """
    Create directory if it doesn't exist
    :param directory: directory to create
    :return: None
    """

    if not os.path.exists(directory):
        os.makedirs(directory)


def generate_file_list(directory, extension):
    """
    Generates a list of files within directory that match extension
    :param directory: directory to search
    :param extension: extension to search for
    :return: a list of paths to files
    """

    counter = 0
    paths = []
    if extension == '.img':
        extension = ('.img')
    elif extension == '.tif':
        extension = ('.tif')

    clear_screen()
    print("Parsing image extens...")

    for root, dirs, files in os.walk(directory):
        for file in files:
            if os.path.splitext(file)[1] in extension:
                path = os.path.join(root, file)
                paths.append(path)
                counter += 1

                if counter % 100 == 0:
                    clear_screen()
                    print("<--- Raster Sectioner --->\n\n")
                    print("Parsed coordinates for {0} images".format(counter))

    return paths


def intersects(image, feature_bounds):
    """
    Checks to see if any of the image lies within the feature
    :param image: list of image boundaries
    :param feature_bounds: list of feature bounds
    :return: Bool
    """

    result = False
    extra_distance = 1000 / 111120  # expand bounding box by this much to ensure overlapping coverage

    feature_ll = (feature_bounds[0] - extra_distance, feature_bounds[1] - extra_distance)
    feature_ur = (feature_bounds[2] + extra_distance, feature_bounds[3] + extra_distance)
    feature_ul = (feature_bounds[0] - extra_distance, feature_bounds[3] + extra_distance)
    feature_lr = (feature_bounds[2] + extra_distance, feature_bounds[1] - extra_distance)

    image_bounds = image[1]

    # Check if any of the image points (corners) lie within, or have a line that crosses, through feature bounding box
    for point in image_bounds:
        if point[0] >= feature_ll[0] and point[0] <= feature_lr[0] and point[1] <= feature_ul[1] and point[1] >= feature_lr[1]:
            result = True

    return result


def write_manifest(path, message):
    """
    Creates a list showing which files went where
    :param path: Path to write file to
    :param message: Message to write
    :return: None
    """

    txtfile = os.path.join(path, 'rs_manifest.txt')
    with open(txtfile, 'a') as f:
        f.write(message)
    f.close()


def clear_screen():
    """
    Clears the terminal window
    :return:
    """

    os.system('cls')


def get_checksum(file):
    """
    Generates SHA1 checksum of file
    :param file: file to check
    :return: stirng
    """

    buffer_size = 256
    with open(file, 'rb') as f:
        while True:
            data = f.read(buffer_size)
            return sha1(data).hexdigest()

def choose_image_type():
    """
    Prompts user to choose image type
    """

    image_type = input("\nImage type:\n"
                       "1):     *.tif\n"
                       "2):     *.img\n"
                       ">>> ")

    if image_type == '1':
        image_type = ".tif"
    elif image_type == '2':
        image_type = ".img"
    else:
        choose_image_type()

    return image_type


def main():
    """ run it all"""

    raster_geometries = []
    raster_functions = Raster()
    shape_functions = Shapes()

    raster_data = input("Raster root directory: ")
    image_type = choose_image_type()
    shape_data = input("Section shapefile: ")
    output_dir = input("Raster output directory: ")

    # Set level to "logging.WARNING" to quieten the log output
    logging.basicConfig(filename=os.path.join(output_dir, 'rs.log'),
                        format='%(asctime)s %(levelname)s %(message)s',
                        level=logging.INFO, filemode='w')

    logging.info("Raster input directory: {0}".format(raster_data))
    logging.info("Raster types: {0}".format(image_type))
    logging.info("Input shapefile: {0}".format(shape_data))
    logging.info("Raster output directory: {0}".format(output_dir))
    logging.info("Searching raster root directory for imagery")

    raster_files = generate_file_list(raster_data, image_type)
    logging.info("Found {0} raster files".format(len(raster_files)))

    for image in raster_files:
        raster_geometries[image] = raster_functions.get_raster_extents(image)
        logging.info("Calculated geometries for {0}".format(image))

    feature_bounds = shape_functions.get_feature_bounds(shape_data)

    for feature in feature_bounds:
        new_directory = os.path.join(os.path.join(output_dir, 'subsections'), feature[0])
        make_directory(new_directory)
        logging.info("Created directory: {0}".format(new_directory))

    for feature in feature_bounds:
        clear_screen()
        print("-> Copying imagery for {0}".format(feature[0]))
        for image, bounds in raster_geometries.items():
            if intersects((image, bounds), feature[1]):
                _output_dir = os.path.join(os.path.join(output_dir, 'subsections'), feature[0])
                output_file = os.path.join(_output_dir, ntpath.basename(image))
                text = "-> Matched raster {0} to subsection {1}".format(image, feature[0])
                logging.info(text)

                try:
                    print(" -> Generating checksum for {0}".format(image))
                    original_checksum = get_checksum(image)
                    if not os.path.exists(output_file) or os.path.exists(output_file) and \
                        original_checksum != get_checksum(output_file): # prevent re-copying and prevent copy errors
                        print(" ->  {0}:    Copying {1} to {2}".format(datetime.datetime.now(), image, output_file))

                        shutil.copyfile(image, output_file)

                        # Ensure output file is identical to input
                        while original_checksum != get_checksum(output_file):
                            logging.warning("File hash mismatch. Attempting copy again. {0}".format(image))
                            print("*** File hash mismatch. Attempting copy again. {0} ***".format(image))
                            shutil.copyfile(image, output_file)

                        write_manifest(image, output_file)

                    else:
                        logging.warning("Raster file {0} already exists. Did not overwrite".format(output_file))

                    # Copy extra data files
                    image_name = ntpath.basename(image)
                    other_file_extensions = ('.ige', '.rrd', '.rde')
                    for extension in other_file_extensions:
                        filename = "{0}{1}".format(os.path.splitext(image_name)[0], extension)
                        path = os.path.dirname(image)
                        extra_file_path = os.path.join(path, filename)
                        if os.path.exists(extra_file_path):
                            extra_checksum = get_checksum(extra_file_path)

                            output_extrafile = os.path.join(_output_dir, filename)
                            print("Copying extrafile {0}".format(extra_file_path))
                            write_manifest(_output_dir, "Copying {0}".format(extra_file_path))

                            if os.path.exists(output_extrafile) and extra_checksum != get_checksum(output_extrafile):
                                if not os.path.exists(output_extrafile):
                                    shutil.copyfile(extra_file_path, output_extrafile)

                                    while extra_checksum != get_checksum(extra_file_path):
                                        logging.warning("File hash mismatch. Attempting copy again. {0]".format(extra_file_path))
                                        print("File hash mismatch. Attempting copy again. {0}".format(extra_file_path))
                                        shutil.copyfile(extra_file_path, output_extrafile)

                                    write_manifest(_output_dir, "\n{0} -----> {1}".format(extra_file_path, output_extrafile))

                                else:
                                    logging.warning("Extra file {0} already exists. Did not overwrite.".format(output_extrafile))

                                logging.info("-> Copied file {0}".format(extra_file_path))

                    text = "-> Copied {0} to {1}".format(image, output_file)
                    logging.info(text)
                    print(text)
                except Exception as e:
                    text = "-> Could not copy {0} to {1}. Returned {2}".format(image, output_file, e)
                    logging.error(text)
                    print(text)

        text = "-> Completed copying imagery for {0}\n".format(feature[0])
        print(text)
        logging.info(text)
    logging.info("-- ** -- COPYING COMPLETE --** --")
    input("-- ** -- COPYING COMPLETE -- ** --")

if __name__ == '__main__':
    main()