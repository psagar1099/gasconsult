#!/usr/bin/env python3
"""
STANDALONE SCRIPT - Run this locally to generate real DALL-E Cormack-Lehane images

REQUIREMENTS:
  pip install openai requests

USAGE:
  export OPENAI_API_KEY="your-api-key-here"
  python GENERATE_DALLE_IMAGES.py

This will generate 4 PNG images in the current directory.
Then upload them to your Render project's static/ directory.
"""

import os
import requests
import time
import json

def main():
    # Check for API key
    api_key = os.getenv("OPENAI_API_KEY")
    if not api_key:
        print("ERROR: Please set OPENAI_API_KEY environment variable")
        print("Example: export OPENAI_API_KEY='sk-...'")
        return
    
    try:
        from openai import OpenAI
        client = OpenAI(api_key=api_key)
    except ImportError:
        print("ERROR: OpenAI library not installed")
        print("Run: pip install openai")
        return
    
    # Medical photography prompts for hyperrealistic endoscopic views
    prompts = {
        1: """Medical endoscopy photograph: Cormack-Lehane Grade 1 laryngoscopy view through videolaryngoscope.
        Hyperrealistic clinical photography showing COMPLETE VISUALIZATION of glottis.
        Sharp, clear view of: bilateral white pearly vocal cords forming V-shape with crisp edges,
        paired coral-pink arytenoid cartilages posteriorly, wide dark triangular glottic opening (rima glottidis),
        surrounding pink-red moist pharyngeal mucosa with visible small blood vessels.
        Professional medical documentation quality with endoscopic LED lighting creating bright central spot,
        shallow depth of field, glistening wet tissue surfaces, anatomically accurate proportions.
        Clinical photograph suitable for anesthesiology textbook. Grade 1 EASY intubation view.""",
        
        2: """Medical endoscopy photograph: Cormack-Lehane Grade 2 laryngoscopy view through videolaryngoscope.
        Hyperrealistic clinical photography showing PARTIAL VISUALIZATION of glottis.
        Clear view of: posterior portions of vocal cords only (posterior commissure visible),
        prominent paired arytenoid cartilages, anterior commissure BLOCKED by curved epiglottis tissue edge,
        partial view of glottic opening, pink-orange epiglottis tissue entering superior frame obscuring anterior view.
        Professional medical documentation with endoscopic lighting, moist glistening tissue,
        visible blood vessel patterns on mucosa. Clinical photograph for medical education.
        Grade 2 MODERATE difficulty intubation view.""",
        
        3: """Medical endoscopy photograph: Cormack-Lehane Grade 3 laryngoscopy view through videolaryngoscope.
        Hyperrealistic clinical photography showing ONLY EPIGLOTTIS VISIBLE, zero glottic structures.
        Dominant feature: large curved epiglottis appearing as thick leaf-shaped pink-coral cartilaginous structure
        with moist glistening surface, tip curving posteriorly completely blocking glottic view.
        NO vocal cords visible anywhere, NO arytenoids visible, NO glottic opening visible.
        Surrounding deep pink-red pharyngeal mucosa with natural tissue texture and blood vessels.
        Professional medical photography with endoscopic lighting showing tissue depth, curvature, moisture.
        Clinical photograph demonstrating anterior larynx anatomy. Grade 3 DIFFICULT intubation view.""",
        
        4: """Medical endoscopy photograph: Cormack-Lehane Grade 4 laryngoscopy view through videolaryngoscope.
        Hyperrealistic clinical photography showing COMPLETE OBSTRUCTION - zero laryngeal landmarks.
        Only visible: homogeneous soft palate tissue mass and posterior pharyngeal wall mucosa
        appearing as pink-red to deep coral soft tissue completely filling field of view with natural mucosal folds.
        NO epiglottis visible, NO glottic structures, NO vocal cords, NO anatomical landmarks identifiable.
        Uniform soft tissue with moist glistening surface, visible small blood vessels, realistic tissue coloration.
        Professional medical photography with endoscopic lighting showing complete visual obstruction.
        Clinical photograph indicating severe anterior larynx or challenging anatomy.
        Grade 4 VERY DIFFICULT intubation view."""
    }
    
    print("=" * 80)
    print("DALL-E CORMACK-LEHANE IMAGE GENERATOR")
    print("=" * 80)
    print("\nGenerating 4 hyperrealistic medical endoscopy images...")
    print("Total time: approximately 2-4 minutes\n")
    print("Cost: ~$0.16 (4 images × $0.04/image for HD quality)")
    print("=" * 80)
    
    urls = {}
    
    # Generate each image
    for grade in range(1, 5):
        print(f"\n[Step {grade}/4] Generating Cormack-Lehane Grade {grade}...")
        print("-" * 80)
        
        try:
            response = client.images.generate(
                model="dall-e-3",
                prompt=prompts[grade],
                size="1024x1024",
                quality="hd",
                n=1,
            )
            
            url = response.data[0].url
            urls[grade] = url
            print(f"✓ Generated successfully")
            print(f"  URL: {url[:70]}...")
            
            # Rate limiting - wait between requests
            if grade < 4:
                print(f"  [Waiting 10 seconds to avoid rate limiting...]")
                time.sleep(10)
                
        except Exception as e:
            print(f"✗ Generation FAILED: {str(e)}")
            urls[grade] = None
    
    # Download images
    print("\n" + "=" * 80)
    print("DOWNLOADING IMAGES")
    print("=" * 80)
    
    output_files = []
    
    for grade, url in urls.items():
        if url:
            print(f"\n[Download {grade}/4] Downloading Grade {grade} PNG...")
            print("-" * 80)
            
            try:
                response = requests.get(url, timeout=60)
                if response.status_code == 200:
                    filename = f"cormack-lehane-grade-{grade}.png"
                    with open(filename, 'wb') as f:
                        f.write(response.content)
                    
                    size_mb = len(response.content) / (1024 * 1024)
                    print(f"✓ Saved: {filename}")
                    print(f"  Size: {size_mb:.2f} MB")
                    output_files.append(filename)
                else:
                    print(f"✗ Download failed (HTTP {response.status_code})")
            except Exception as e:
                print(f"✗ Download error: {str(e)}")
        else:
            print(f"\n[Download {grade}/4] Skipping Grade {grade} (generation failed)")
    
    # Summary
    print("\n" + "=" * 80)
    print("GENERATION COMPLETE")
    print("=" * 80)
    print(f"\n✓ Successfully created {len(output_files)}/4 images")
    
    if output_files:
        print("\nGenerated files:")
        for filename in output_files:
            print(f"  • {filename}")
        
        print("\n" + "-" * 80)
        print("NEXT STEPS:")
        print("-" * 80)
        print("1. Upload these 4 PNG files to your server's static/ directory")
        print("2. Or commit them to git: git add static/*.png && git commit && git push")
        print("3. Deploy to Render - images will load instantly!")
        print("=" * 80)
    else:
        print("\n✗ No images were successfully generated")
        print("Please check your API key and try again")
    
if __name__ == "__main__":
    main()
