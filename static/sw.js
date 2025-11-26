// Service Worker for gasconsult.ai
// Enables offline functionality for Quick Dose page

const CACHE_NAME = 'gasconsult-v1';
const OFFLINE_CACHE = [
    '/quick-dose',
    '/static/favicon.svg',
    '/static/manifest.json',
    'https://fonts.googleapis.com/css2?family=Sora:wght@300;400;500;600;700&family=Inter:wght@400;500;600&family=JetBrains+Mono:wght@400;500;600&display=swap'
];

// Install event - cache offline pages
self.addEventListener('install', (event) => {
    event.waitUntil(
        caches.open(CACHE_NAME)
            .then((cache) => {
                console.log('Opened cache');
                return cache.addAll(OFFLINE_CACHE);
            })
            .then(() => self.skipWaiting())
    );
});

// Activate event - clean up old caches
self.addEventListener('activate', (event) => {
    event.waitUntil(
        caches.keys().then((cacheNames) => {
            return Promise.all(
                cacheNames.map((cacheName) => {
                    if (cacheName !== CACHE_NAME) {
                        console.log('Deleting old cache:', cacheName);
                        return caches.delete(cacheName);
                    }
                })
            );
        }).then(() => self.clients.claim())
    );
});

// Fetch event - serve from cache when offline
self.addEventListener('fetch', (event) => {
    const url = new URL(event.request.url);

    // Only handle GET requests
    if (event.request.method !== 'GET') {
        return;
    }

    // Network-first strategy for main site
    // Cache-first strategy for Quick Dose page (offline functionality)
    if (url.pathname === '/quick-dose') {
        event.respondWith(
            caches.match(event.request)
                .then((response) => {
                    // Return cached version if available
                    if (response) {
                        // Update cache in background
                        fetch(event.request).then((freshResponse) => {
                            if (freshResponse && freshResponse.status === 200) {
                                caches.open(CACHE_NAME).then((cache) => {
                                    cache.put(event.request, freshResponse.clone());
                                });
                            }
                        });
                        return response;
                    }

                    // If not in cache, fetch from network
                    return fetch(event.request).then((response) => {
                        // Cache the response for future offline use
                        if (response && response.status === 200) {
                            const responseToCache = response.clone();
                            caches.open(CACHE_NAME).then((cache) => {
                                cache.put(event.request, responseToCache);
                            });
                        }
                        return response;
                    });
                })
                .catch(() => {
                    // Return a custom offline page if needed
                    return new Response('Offline - Quick Dose reference will be available once cached', {
                        status: 503,
                        statusText: 'Service Unavailable',
                        headers: new Headers({
                            'Content-Type': 'text/plain'
                        })
                    });
                })
        );
    } else {
        // Network-first for other pages
        event.respondWith(
            fetch(event.request)
                .then((response) => {
                    // Clone the response
                    const responseToCache = response.clone();

                    // Cache successful responses
                    if (response && response.status === 200) {
                        caches.open(CACHE_NAME).then((cache) => {
                            cache.put(event.request, responseToCache);
                        });
                    }

                    return response;
                })
                .catch(() => {
                    // If network fails, try cache
                    return caches.match(event.request);
                })
        );
    }
});
