// Node.js script to generate a manifest of available scenes
// Run with: node generate_manifest.js

const fs = require('fs');
const path = require('path');

const dataDir = './pgsr_data';
const outputFile = './scenes_manifest.json';

function parseFilename(filename) {
    // Parse filenames like "gs_aomlion_lion1.splat" or "mesh_dtu_dtu-scan105.ply"
    const match = filename.match(/^(gs|mesh)_([^_]+)_(.+)\.(splat|ply)$/);
    if (!match) return null;
    
    return {
        type: match[1], // 'gs' or 'mesh'
        dataset: match[2],
        sceneName: match[3].replace(/\.(splat|ply)$/, ''),
        extension: match[4]
    };
}

function generateManifest() {
    if (!fs.existsSync(dataDir)) {
        console.error(`Directory ${dataDir} does not exist`);
        process.exit(1);
    }

    const files = fs.readdirSync(dataDir);
    const scenes = {};

    files.forEach(file => {
        const parsed = parseFilename(file);
        if (!parsed) return;

        const sceneKey = `${parsed.dataset}_${parsed.sceneName}`;
        
        if (!scenes[sceneKey]) {
            scenes[sceneKey] = {
                dataset: parsed.dataset,
                sceneName: parsed.sceneName,
                displayName: `${parsed.dataset} - ${parsed.sceneName}`,
                splat: null,
                mesh: null
            };
        }

        if (parsed.type === 'gs') {
            if (parsed.extension === 'splat') {
                scenes[sceneKey].splat = file;
            }
        } else if (parsed.type === 'mesh') {
            scenes[sceneKey].mesh = file;
        }
    });

    // Filter to only include scenes with both .splat and mesh.ply
    const validScenes = Object.entries(scenes)
        .filter(([key, scene]) => scene.splat && scene.mesh)
        .reduce((acc, [key, scene]) => {
            acc[key] = scene;
            return acc;
        }, {});

    const manifest = {
        generatedAt: new Date().toISOString(),
        dataDir: 'pgsr_data',
        scenes: validScenes
    };

    fs.writeFileSync(outputFile, JSON.stringify(manifest, null, 2));
    console.log(`Manifest generated with ${Object.keys(validScenes).length} scenes`);
    console.log(`Written to ${outputFile}`);
}

generateManifest();
