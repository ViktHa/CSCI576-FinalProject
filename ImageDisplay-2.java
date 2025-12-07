import java.awt.*;
import java.awt.image.*;
import java.io.*;
import java.util.*;
import java.util.List;
import javax.swing.*;

public class ImageDisplay {
    JFrame frame;
    JLabel lbIm1;
    BufferedImage imgOne;
    int width;
    int height;

    String imagePath;

    // Detected puzzle pieces
    List<PuzzlePiece> pieces = new ArrayList<>();

    // Class to represent a detected puzzle piece
    static class PuzzlePiece {
        BufferedImage image;   // final, axis-aligned, square tile
        Rectangle bounds;      // original bounding box in source image
        int id;
        double rotationAngle;  // angle (in degrees) applied to make it axis-aligned

        PuzzlePiece(BufferedImage image, Rectangle bounds, int id) {
            this.image = image;
            this.bounds = bounds;
            this.id = id;
            this.rotationAngle = 0;
        }
    }

    public void readImageRGB(int width, int height, String imgPath, BufferedImage img) {
        try {
            int frameLength = width * height * 3;
            File file = new File(imgPath);
            RandomAccessFile raf = new RandomAccessFile(file, "r");
            raf.seek(0);
            byte[] bytes = new byte[frameLength];
            raf.read(bytes);

            int ind = 0;
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    byte r = bytes[ind];
                    byte g = bytes[ind + height * width];
                    byte b = bytes[ind + height * width * 2];
                    int pix = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
                    img.setRGB(x, y, pix);
                    ind++;
                }
            }
            raf.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void parseArgs(String[] args) {
        imagePath = args[0];
        width = Integer.parseInt(args[1]);
        height = Integer.parseInt(args[2]);
    }

    /**
     * Detects puzzle pieces using connected component analysis
     * Assumes background is relatively uniform (e.g., black)
     */
    public void detectPieces() {
        // Create binary mask: true = foreground (piece), false = background
        boolean[][] mask = createForegroundMask();

        // Find connected components
        boolean[][] visited = new boolean[height][width];
        int pieceId = 0;

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                if (mask[y][x] && !visited[y][x]) {
                    // Found a new piece, extract it
                    PuzzlePiece piece = extractPiece(x, y, mask, visited, pieceId);
                    if (piece != null && piece.bounds.width > 10 && piece.bounds.height > 10) {
                        pieces.add(piece);
                        pieceId++;
                    }
                }
            }
        }

        System.out.println("Detected " + pieces.size() + " puzzle pieces");

        // Print rotation information
        for (PuzzlePiece piece : pieces) {
            System.out.println("Piece " + piece.id + " - Rotation used: " + piece.rotationAngle + " degrees");
        }
    }

    /**
     * Creates a binary mask separating foreground (pieces) from background
     */
    private boolean[][] createForegroundMask() {
        boolean[][] mask = new boolean[height][width];

        // Sample background color from a corner
        int bgColor = estimateBackgroundColor();
        int threshold = 30; // Adjust if needed

        int bgR = (bgColor >> 16) & 0xFF;
        int bgG = (bgColor >> 8) & 0xFF;
        int bgB = bgColor & 0xFF;

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int rgb = imgOne.getRGB(x, y);
                int r = (rgb >> 16) & 0xFF;
                int g = (rgb >> 8) & 0xFF;
                int b = rgb & 0xFF;

                // Check if pixel is significantly different from background
                int diff = Math.abs(r - bgR) + Math.abs(g - bgG) + Math.abs(b - bgB);
                mask[y][x] = diff > threshold;
            }
        }

        return mask;
    }

    /**
     * Estimates background color by sampling corners
     */
    private int estimateBackgroundColor() {
        int[] corners = new int[4];
        corners[0] = imgOne.getRGB(0, 0);
        corners[1] = imgOne.getRGB(width - 1, 0);
        corners[2] = imgOne.getRGB(0, height - 1);
        corners[3] = imgOne.getRGB(width - 1, height - 1);

        // Simple: just use the first corner (your backgrounds are uniform)
        return corners[0];
    }

    /**
     * Extracts a single piece using flood fill to find connected region
     */
    private PuzzlePiece extractPiece(int startX, int startY, boolean[][] mask, boolean[][] visited, int id) {
        Queue<Point> queue = new LinkedList<>();
        List<Point> piecePixels = new ArrayList<>();

        queue.add(new Point(startX, startY));
        visited[startY][startX] = true;

        int minX = startX, maxX = startX;
        int minY = startY, maxY = startY;

        // Flood fill to find all pixels in this piece
        while (!queue.isEmpty()) {
            Point p = queue.poll();
            piecePixels.add(p);

            minX = Math.min(minX, p.x);
            maxX = Math.max(maxX, p.x);
            minY = Math.min(minY, p.y);
            maxY = Math.max(maxY, p.y);

            // Check 8-connected neighbors
            for (int dy = -1; dy <= 1; dy++) {
                for (int dx = -1; dx <= 1; dx++) {
                    if (dx == 0 && dy == 0) continue;

                    int nx = p.x + dx;
                    int ny = p.y + dy;

                    if (nx >= 0 && nx < width && ny >= 0 && ny < height &&
                            mask[ny][nx] && !visited[ny][nx]) {
                        visited[ny][nx] = true;
                        queue.add(new Point(nx, ny));
                    }
                }
            }
        }

        // Local coordinates relative to bounding box
        Rectangle bounds = new Rectangle(minX, minY, maxX - minX + 1, maxY - minY + 1);
        List<Point> localPixels = new ArrayList<>();

        for (Point p : piecePixels) {
            int lx = p.x - bounds.x;
            int ly = p.y - bounds.y;
            localPixels.add(new Point(lx, ly));
        }

        // Detect rotation angle by brute-force search (min-area bounding box)
        double rotationAngle = detectRotation(localPixels);

        // Extract the piece and rotate it to axis-aligned
        BufferedImage pieceImage = extractAndRotatePiece(localPixels, bounds, rotationAngle);

        PuzzlePiece piece = new PuzzlePiece(pieceImage, bounds, id);
        piece.rotationAngle = rotationAngle;

        return piece;
    }

    /**
     * Brute-force rotation detection:
     * For angles in [-89, 89], rotate all local pixels and find the
     * bounding box. The angle that minimizes area corresponds to an
     * axis-aligned orientation of the square piece.
     *
     * The returned angle is what we should apply to the subimage.
     */
    private double detectRotation(List<Point> pixels) {
        int n = pixels.size();
        double[] xs = new double[n];
        double[] ys = new double[n];

        for (int i = 0; i < n; i++) {
            xs[i] = pixels.get(i).x;
            ys[i] = pixels.get(i).y;
        }

        double bestAngle = 0.0;
        double bestArea = Double.MAX_VALUE;

        // 1 degree resolution is enough; you can reduce step if desired.
        for (double angle = -89.0; angle <= 89.0; angle += 1.0) {
            double rad = Math.toRadians(angle);
            double cos = Math.cos(rad);
            double sin = Math.sin(rad);

            double minX = Double.POSITIVE_INFINITY;
            double maxX = Double.NEGATIVE_INFINITY;
            double minY = Double.POSITIVE_INFINITY;
            double maxY = Double.NEGATIVE_INFINITY;

            for (int i = 0; i < n; i++) {
                double rx = xs[i] * cos - ys[i] * sin;
                double ry = xs[i] * sin + ys[i] * cos;

                if (rx < minX) minX = rx;
                if (rx > maxX) maxX = rx;
                if (ry < minY) minY = ry;
                if (ry > maxY) maxY = ry;
            }

            double w = maxX - minX;
            double h = maxY - minY;
            double area = w * h;

            if (area < bestArea) {
                bestArea = area;
                bestAngle = angle;
            }
        }

        return bestAngle;
    }

    /**
     * Extracts piece and rotates it to correct orientation.
     * pixels: LOCAL coordinates relative to 0..bounds.width/height
     */
    private BufferedImage extractAndRotatePiece(
            List<Point> pixels, Rectangle bounds, double angleDegrees) {

        // Extract local subimage (ARGB, transparent background)
        BufferedImage original = new BufferedImage(bounds.width, bounds.height, BufferedImage.TYPE_INT_ARGB);

        for (Point p : pixels) {
            int lx = p.x;
            int ly = p.y;
            int sx = bounds.x + lx;
            int sy = bounds.y + ly;
            original.setRGB(lx, ly, imgOne.getRGB(sx, sy));
        }

        // If rotation is ~0, just crop to content and pad to square
        if (Math.abs(angleDegrees) < 1e-3) {
            BufferedImage cropped0 = cropToContent(original);
            return padToSquare(cropped0, Color.BLACK);
        }

        double angle = Math.toRadians(angleDegrees);

        // Compute rotated image bounds
        double sin = Math.abs(Math.sin(angle));
        double cos = Math.abs(Math.cos(angle));
        int newW = (int) Math.ceil(bounds.width * cos + bounds.height * sin);
        int newH = (int) Math.ceil(bounds.height * cos + bounds.width * sin);

        BufferedImage rotated = new BufferedImage(newW, newH, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g2d = rotated.createGraphics();
        g2d.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);

        // Clear with transparency
        g2d.setComposite(AlphaComposite.Clear);
        g2d.fillRect(0, 0, newW, newH);
        g2d.setComposite(AlphaComposite.Src);

        // Rotate around the center of the original bounding box
        g2d.translate(newW / 2.0, newH / 2.0);
        g2d.rotate(angle);
        g2d.translate(-bounds.width / 2.0, -bounds.height / 2.0);
        g2d.drawImage(original, 0, 0, null);
        g2d.dispose();

        // Crop away all transparent triangles
        BufferedImage cropped = cropToContent(rotated);

        // Pad to centered square with black background (RGB)
        return padToSquare(cropped, Color.BLACK);
    }

    /**
     * Crops an ARGB image to its non-transparent content.
     */
    private BufferedImage cropToContent(BufferedImage img) {
        int w = img.getWidth();
        int h = img.getHeight();

        int minX = w, minY = h, maxX = 0, maxY = 0;
        boolean found = false;

        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                int a = (img.getRGB(x, y) >> 24) & 0xFF;
                if (a > 10) { // non-transparent
                    found = true;
                    if (x < minX) minX = x;
                    if (x > maxX) maxX = x;
                    if (y < minY) minY = y;
                    if (y > maxY) maxY = y;
                }
            }
        }

        if (!found) {
            return img;
        }

        return img.getSubimage(minX, minY, maxX - minX + 1, maxY - minY + 1);
    }

    /**
     * Pads an image to a square of max(width,height) with a solid background,
     * centering the content.
     */
    private BufferedImage padToSquare(BufferedImage img, Color bg) {
        int w = img.getWidth();
        int h = img.getHeight();
        int size = Math.max(w, h);

        BufferedImage sq = new BufferedImage(size, size, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = sq.createGraphics();
        g.setColor(bg);
        g.fillRect(0, 0, size, size);

        int dx = (size - w) / 2;
        int dy = (size - h) / 2;
        g.drawImage(img, dx, dy, null);
        g.dispose();

        return sq;
    }

    /**
     * Shows the rotated pieces composited back into a blank canvas
     * (optional, but handy for visually checking orientation).
     */
    public void showPiecesRotatedOnCanvas() {
        BufferedImage display = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = display.createGraphics();
        g.setColor(Color.BLACK);
        g.fillRect(0, 0, width, height);

        for (PuzzlePiece piece : pieces) {
            BufferedImage rotated = piece.image;
            int px = piece.bounds.x;
            int py = piece.bounds.y;

            // center the rotated piece at its original bounding-box center
            int drawX = px + (piece.bounds.width - rotated.getWidth()) / 2;
            int drawY = py + (piece.bounds.height - rotated.getHeight()) / 2;

            g.drawImage(rotated, drawX, drawY, null);
        }

        g.dispose();
        lbIm1 = new JLabel(new ImageIcon(display));
    }

    /**
     * Displays all extracted pieces in a grid.
     * Each piece.image is already axis-aligned and square with no triangles.
     */
    public void showExtractedPieces() {
        if (pieces.isEmpty()) {
            System.out.println("No pieces detected!");
            return;
        }

        JFrame pieceFrame = new JFrame("Extracted Pieces");
        pieceFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        JPanel panel = new JPanel();
        panel.setLayout(new GridLayout(0, 4, 10, 10)); // 4 columns

        for (PuzzlePiece piece : pieces) {
            JLabel label = new JLabel(new ImageIcon(piece.image));
            label.setBorder(BorderFactory.createTitledBorder("Piece " + piece.id));
            panel.add(label);
        }

        JScrollPane scrollPane = new JScrollPane(panel);
        pieceFrame.add(scrollPane);
        pieceFrame.setSize(800, 600);
        pieceFrame.setVisible(true);
    }

    // Edge types for matching
    enum EdgeType {
        TOP, BOTTOM, LEFT, RIGHT
    }

    // Represents an edge match between two pieces
    static class EdgeMatch implements Comparable<EdgeMatch> {
        int piece1Id;
        EdgeType edge1;
        int piece2Id;
        EdgeType edge2;
        double score; // Lower is better

        EdgeMatch(int p1, EdgeType e1, int p2, EdgeType e2, double score) {
            this.piece1Id = p1;
            this.edge1 = e1;
            this.piece2Id = p2;
            this.edge2 = e2;
            this.score = score;
        }

        @Override
        public int compareTo(EdgeMatch other) {
            return Double.compare(this.score, other.score);
        }

        @Override
        public String toString() {
            return String.format("Piece %d (%s) <-> Piece %d (%s): %.2f",
                    piece1Id, edge1, piece2Id, edge2, score);
        }
    }

    /**
     * Extracts pixels along a specific edge of a piece
     */
    private int[][] extractEdgePixels(PuzzlePiece piece, EdgeType edge, int stripWidth) {
        BufferedImage img = piece.image;  // already axis-aligned & square
        int w = img.getWidth();
        int h = img.getHeight();

        int[][] edgePixels = null;

        switch (edge) {
            case TOP:
                int topHeight = Math.min(stripWidth, h);
                edgePixels = new int[w][topHeight];
                for (int x = 0; x < w; x++) {
                    for (int y = 0; y < topHeight; y++) {
                        edgePixels[x][y] = img.getRGB(x, y);
                    }
                }
                break;

            case BOTTOM:
                int bottomHeight = Math.min(stripWidth, h);
                edgePixels = new int[w][bottomHeight];
                for (int x = 0; x < w; x++) {
                    for (int y = 0; y < bottomHeight; y++) {
                        edgePixels[x][y] = img.getRGB(x, h - bottomHeight + y);
                    }
                }
                break;

            case LEFT:
                int leftWidth = Math.min(stripWidth, w);
                edgePixels = new int[h][leftWidth];
                for (int y = 0; y < h; y++) {
                    for (int x = 0; x < leftWidth; x++) {
                        edgePixels[y][x] = img.getRGB(x, y);
                    }
                }
                break;

            case RIGHT:
                int rightWidth = Math.min(stripWidth, w);
                edgePixels = new int[h][rightWidth];
                for (int y = 0; y < h; y++) {
                    for (int x = 0; x < rightWidth; x++) {
                        edgePixels[y][x] = img.getRGB(w - rightWidth + x, y);
                    }
                }
                break;
        }

        return edgePixels;
    }

    /**
     * Compares two edges and computes compatibility score using pixel comparison
     */
    private double compareEdges(int[][] edge1, int[][] edge2) {
        if (edge1.length != edge2.length) {
            return Double.MAX_VALUE; // Incompatible sizes
        }

        int length = edge1.length;
        int depth = Math.min(edge1[0].length, edge2[0].length);

        double totalError = 0;
        int pixelCount = 0;

        for (int i = 0; i < length; i++) {
            for (int d = 0; d < depth; d++) {
                int rgb1 = edge1[i][d];
                int rgb2 = edge2[i][d];

                int r1 = (rgb1 >> 16) & 0xFF;
                int g1 = (rgb1 >> 8) & 0xFF;
                int b1 = rgb1 & 0xFF;

                int r2 = (rgb2 >> 16) & 0xFF;
                int g2 = (rgb2 >> 8) & 0xFF;
                int b2 = rgb2 & 0xFF;

                totalError += Math.abs(r1 - r2) + Math.abs(g1 - g2) + Math.abs(b1 - b2);
                pixelCount++;
            }
        }

        return totalError / pixelCount;
    }

    private EdgeType getComplementaryEdge(EdgeType edge) {
        switch (edge) {
            case TOP:
                return EdgeType.BOTTOM;
            case BOTTOM:
                return EdgeType.TOP;
            case LEFT:
                return EdgeType.RIGHT;
            case RIGHT:
                return EdgeType.LEFT;
            default:
                return null;
        }
    }

    /**
     * Performs edge matching for all piece pairs
     */
    public List<EdgeMatch> performEdgeMatching(int stripWidth) {
        List<EdgeMatch> matches = new ArrayList<>();

        System.out.println("Starting edge matching...");

        for (int i = 0; i < pieces.size(); i++) {
            for (int j = i + 1; j < pieces.size(); j++) {
                PuzzlePiece piece1 = pieces.get(i);
                PuzzlePiece piece2 = pieces.get(j);

                for (EdgeType edge1 : EdgeType.values()) {
                    EdgeType edge2 = getComplementaryEdge(edge1);

                    int[][] pixels1 = extractEdgePixels(piece1, edge1, stripWidth);
                    int[][] pixels2 = extractEdgePixels(piece2, edge2, stripWidth);

                    double score = compareEdges(pixels1, pixels2);
                    if (score != Double.MAX_VALUE) {
                        matches.add(new EdgeMatch(piece1.id, edge1, piece2.id, edge2, score));
                    }
                }
            }
        }

        Collections.sort(matches);

        System.out.println("Found " + matches.size() + " potential matches");
        System.out.println("\nTop 10 matches:");
        for (int i = 0; i < Math.min(10, matches.size()); i++) {
            System.out.println(matches.get(i));
        }

        return matches;
    }

    public void showIms() {
        frame = new JFrame();
        GridBagLayout gLayout = new GridBagLayout();
        frame.getContentPane().setLayout(gLayout);

        GridBagConstraints c = new GridBagConstraints();
        c.gridx = 0;
        c.gridy = 0;

        frame.getContentPane().add(lbIm1, c);

        frame.pack();
        frame.setVisible(true);
    }

    public static void main(String[] args) {
        ImageDisplay iq = new ImageDisplay();
        iq.parseArgs(args);

        iq.imgOne = new BufferedImage(iq.width, iq.height, BufferedImage.TYPE_INT_RGB);
        iq.readImageRGB(iq.width, iq.height, iq.imagePath, iq.imgOne);

        // Detect and extract pieces (with rotation normalization)
        iq.detectPieces();

        // Show rotated pieces composited on a canvas (optional)
        iq.showPiecesRotatedOnCanvas();
        iq.showIms();

        // Show extracted axis-aligned pieces in a separate window
        iq.showExtractedPieces();

        // Perform edge matching on normalized tiles
        int stripWidth = 5;
        iq.performEdgeMatching(stripWidth);
    }
}
