# DALL-E Cormack-Lehane Image Generation

## Overview

The Difficult Airway Assessment feature uses **real DALL-E generated medical endoscopy photographs** showing each Cormack-Lehane grade. These are hyperrealistic AI-generated images that look like actual laryngoscopy views.

## Current Status

The code is configured to use PNG images at:
```
/static/cormack-lehane-grade-1.png  (Grade 1: Full glottic view - EASY)
/static/cormack-lehane-grade-2.png  (Grade 2: Partial view - MODERATE)
/static/cormack-lehane-grade-3.png  (Grade 3: Epiglottis only - DIFFICULT)
/static/cormack-lehane-grade-4.png  (Grade 4: No structures - VERY DIFFICULT)
```

## How to Generate Real DALL-E Images

### Option 1: Run the Python Script (Recommended)

1. **Install requirements:**
   ```bash
   pip install openai requests
   ```

2. **Set your OpenAI API key:**
   ```bash
   export OPENAI_API_KEY="sk-your-key-here"
   ```

3. **Run the generation script:**
   ```bash
   python GENERATE_DALLE_IMAGES.py
   ```

4. **Wait 2-4 minutes** for all 4 images to generate

5. **Upload the PNG files:**
   ```bash
   # Copy to static directory
   cp cormack-lehane-grade-*.png static/
   
   # Or commit and push
   git add static/cormack-lehane-grade-*.png
   git commit -m "Add DALL-E generated Cormack-Lehane images"
   git push
   ```

### Option 2: Manual Generation via OpenAI Dashboard

1. Go to https://platform.openai.com/playground/images
2. Use the prompts from `GENERATE_DALLE_IMAGES.py` (lines 28-72)
3. Select model: `DALL-E 3`, size: `1024x1024`, quality: `HD`
4. Generate each image and download
5. Rename files to match the expected names
6. Upload to `static/` directory

### Option 3: Use Temporary SVG Placeholders

If you want to test without DALL-E images, temporary SVG illustrations are available:
```bash
# Already in static/ directory
static/cormack-lehane-grade-1.svg
static/cormack-lehane-grade-2.svg
static/cormack-lehane-grade-3.svg
static/cormack-lehane-grade-4.svg
```

**Note:** Update `app.py` line 28526 to use `.svg` instead of `.png` if using placeholders.

## Cost

- **DALL-E 3 HD quality:** $0.04 per image
- **Total for 4 images:** $0.16 (one-time cost)

## Technical Details

### Image Specifications
- **Model:** DALL-E 3
- **Quality:** HD
- **Size:** 1024×1024 pixels
- **Format:** PNG
- **File size:** ~1-3 MB per image
- **Style:** Hyperrealistic medical endoscopy photography

### Prompt Design
The prompts are designed to generate:
- Clinical photography quality (not illustrations/cartoons)
- Accurate anatomical landmarks for each grade
- Realistic endoscopic lighting and tissue appearance
- Moist glistening mucosa surfaces
- Proper tissue colors (pink-red pharynx, white vocal cords, coral epiglottis)
- Professional medical documentation style

### Why Static Images?

1. **Performance:** Instant loading (no 30+ second DALL-E generation time)
2. **Reliability:** No worker timeouts or API failures
3. **Consistency:** Same high-quality image every time
4. **Cost:** One-time $0.16 vs $0.04 per assessment
5. **Standardization:** Cormack-Lehane grades are standardized anyway

## Troubleshooting

**Q: Images not appearing on website?**
- Verify PNG files exist in `/static/` directory
- Check file names match exactly: `cormack-lehane-grade-{1-4}.png`
- Clear browser cache and hard refresh (Cmd+Shift+R or Ctrl+Shift+F5)

**Q: DALL-E generation failed?**
- Check your OpenAI API key is valid and has credits
- Ensure you have DALL-E API access enabled
- Wait a few seconds and try again (rate limiting)
- Check OpenAI status page for outages

**Q: Want to regenerate with better prompts?**
- Edit the prompts in `GENERATE_DALLE_IMAGES.py`
- Re-run the script
- Replace the old PNG files

## Future Enhancements

Potential improvements:
- [ ] Generate images with different lighting conditions
- [ ] Create variations for edge cases (blood, secretions, etc.)
- [ ] Add 4K resolution versions
- [ ] Implement lazy loading for faster page performance
- [ ] Add image compression to reduce file sizes

---

**Last Updated:** December 21, 2025
